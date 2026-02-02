import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

class LaplaceBEM2D:
    def __init__(self, n_elements_per_side=20):
        """
        Initializes the BEM solver for a unit square [0,1]x[0,1].
        Uses Constant Elements to avoid corner singularities.
        """
        self.N_side = n_elements_per_side
        self.N_total = 4 * n_elements_per_side
        
        # Mesh geometry
        self.nodes = None      # Collocation points (midpoints)
        self.endpoints = None  # (x1, x2, y1, y2) for each element
        self.normals = None    # Normal vectors (nx, ny)
        self.lengths = None    # Element lengths
        
        # Solution storage
        self.u_boundary = None # Known Dirichlet values
        self.q_boundary = None # Computed Neumann fluxes
        
        self._generate_mesh()

    def _generate_mesh(self):
        """Generates mesh counter-clockwise starting from (0,0)."""
        corners = [(0,0), (1,0), (1,1), (0,1)] # Bottom, Right, Top, Left
        
        x_ends, y_ends = [], []
        
        for i in range(4):
            start, end = corners[i], corners[(i+1)%4]
            xs = np.linspace(start[0], end[0], self.N_side + 1)
            ys = np.linspace(start[1], end[1], self.N_side + 1)
            
            for k in range(self.N_side):
                x_ends.append([xs[k], xs[k+1]])
                y_ends.append([ys[k], ys[k+1]])
                
        self.endpoints = np.column_stack((np.array(x_ends), np.array(y_ends)))
        
        # Compute geometric properties (nodes, lengths, normals)
        self.nodes = np.zeros((self.N_total, 2))
        self.normals = np.zeros((self.N_total, 2))
        self.lengths = np.zeros(self.N_total)
        
        for i in range(self.N_total):
            x1, x2 = self.endpoints[i, 0], self.endpoints[i, 1]
            y1, y2 = self.endpoints[i, 2], self.endpoints[i, 3]
            
            # Midpoint (Node)
            self.nodes[i] = [(x1+x2)/2, (y1+y2)/2]
            
            # Length and Normal
            dx, dy = x2 - x1, y2 - y1
            self.lengths[i] = np.sqrt(dx**2 + dy**2)
            self.normals[i] = [dy/self.lengths[i], -dx/self.lengths[i]]

    def _fundamental_solution(self, r, n_y, r_vec):
        """
        Computes Green's function u* and its derivative q*.
        u* = (1/2pi) * ln(1/r)
        q* = -(1/2pi) * (r . n) / r^2
        """
        const = 1.0 / (2 * np.pi)
        if r < 1e-15: return 0, 0 # Handled analytically in diagonal terms
        
        u_star = const * np.log(1.0/r)
        dr_dn = np.dot(r_vec, n_y)
        q_star = -const * dr_dn / (r**2)
        return u_star, q_star

    def solve(self, dirichlet_func):
        """
        Solves the BEM system H u = G q.
        Args:
            dirichlet_func: A python function f(x,y) that returns the u value.
        """
        # 1. Apply Functional BCs to all boundary nodes
        self.u_boundary = np.array([dirichlet_func(n[0], n[1]) for n in self.nodes])
        
        # 2. Assemble Matrices H and G
        H = np.zeros((self.N_total, self.N_total))
        G = np.zeros((self.N_total, self.N_total))
        
        # Gauss quadrature constants (4-point rule)
        gauss_pt = np.array([-0.861136, -0.339981, 0.339981, 0.861136])
        gauss_wt = np.array([0.347855, 0.652145, 0.652145, 0.347855])
        
        print(f"Assembling matrices ({self.N_total}x{self.N_total})...")
        
        for i in range(self.N_total): # Source Node i
            source = self.nodes[i]
            for j in range(self.N_total): # Field Element j
                if i == j:
                    # Singular (Diagonal) terms
                    # H_ii: For constant elements on smooth boundary, integral is 0. 
                    #       Remaining term is c=0.5 (free term).
                    H[i, j] = 0.5 
                    
                    # G_ii: Analytical integral of ln(1/r)
                    L = self.lengths[j]
                    G[i, j] = (L / (2*np.pi)) * (np.log(2.0/L) + 1.0)
                else:
                    # Regular (Off-diagonal) terms -> Gauss Integration
                    x1, x2 = self.endpoints[j, 0], self.endpoints[j, 1]
                    y1, y2 = self.endpoints[j, 2], self.endpoints[j, 3]
                    L = self.lengths[j]
                    n_j = self.normals[j]
                    
                    int_u, int_q = 0.0, 0.0
                    for k in range(4):
                        xi = gauss_pt[k]
                        # Interpolate to global coordinates
                        x_g = 0.5*(x1+x2) + 0.5*(x2-x1)*xi
                        y_g = 0.5*(y1+y2) + 0.5*(y2-y1)*xi
                        
                        r_vec = np.array([x_g, y_g]) - source
                        r = np.linalg.norm(r_vec)
                        
                        us, qs = self._fundamental_solution(r, n_j, r_vec)
                        int_u += us * gauss_wt[k]
                        int_q += qs * gauss_wt[k]
                        
                    G[i, j] = int_u * (L/2.0)
                    H[i, j] = int_q * (L/2.0)
        
        # 3. Solve Linear System G * q = H * u
        # Since we have full Dirichlet, u is known, q is unknown.
        # RHS vector b = H * u
        b = np.dot(H, self.u_boundary)
        
        # Solve for q
        self.q_boundary = np.linalg.solve(G, b)
        np.set_printoptions(linewidth=120)
        print(G)
        print(b)
        print(self.q_boundary)

        print("System solved.")

    def evaluate_internal(self, X, Y):
        """Computes potential u at internal grid points X, Y."""
        points = np.column_stack((X.flatten(), Y.flatten()))
        u_internal = np.zeros(len(points))
        
        gauss_pt = np.array([-0.861136, -0.339981, 0.339981, 0.861136])
        gauss_wt = np.array([0.347855, 0.652145, 0.652145, 0.347855])
        
        for k, pt in enumerate(points):
            val = 0.0
            for j in range(self.N_total):
                x1, x2 = self.endpoints[j, 0], self.endpoints[j, 1]
                y1, y2 = self.endpoints[j, 2], self.endpoints[j, 3]
                L = self.lengths[j]
                n_j = self.normals[j]
                u_j = self.u_boundary[j]
                q_j = self.q_boundary[j]
                
                int_G, int_H = 0.0, 0.0
                for gp in range(4):
                    xi = gauss_pt[gp]
                    x_g = 0.5*(x1+x2) + 0.5*(x2-x1)*xi
                    y_g = 0.5*(y1+y2) + 0.5*(y2-y1)*xi
                    
                    r_vec = np.array([x_g, y_g]) - pt
                    r = np.linalg.norm(r_vec)
                    
                    us, qs = self._fundamental_solution(r, n_j, r_vec)
                    int_G += us * gauss_wt[gp]
                    int_H += qs * gauss_wt[gp]
                
                # Internal Formula: u(x) = Integral(u* q) - Integral(q* u)
                val += (int_G * (L/2.0) * q_j) - (int_H * (L/2.0) * u_j)
            u_internal[k] = val
            
        return u_internal.reshape(X.shape)

# --- Define Analytic Solution ---
# We use a known harmonic function to set BCs and verify error.
# u(x,y) = e^x * cos(y) is harmonic (Laplacian = 0).
def harmonic_function(x, y):
    return (x - 0.5) ** 2 - (y - 0.5) ** 2

# --- Main Execution ---

# 1. Setup Solver
bem = LaplaceBEM2D(n_elements_per_side=2)

# 2. Solve using the function
# We pass the harmonic function itself to serve as the Dirichlet BC generator
bem.solve(harmonic_function)

# 3. Create Internal Grid for Plotting
# Note: avoid boundary edges (0.0 and 1.0) strictly to avoid singularities in plotting
grid_n = 20
x_space = np.linspace(0.05, 0.95, grid_n)
y_space = np.linspace(0.05, 0.95, grid_n)
X, Y = np.meshgrid(x_space, y_space)

print("Computing internal solution...")
U_bem = bem.evaluate_internal(X, Y)

# 4. Compute Exact Solution and Error
U_exact = harmonic_function(X, Y)
U_error = np.abs(U_bem - U_exact)

# --- Visualization Control ---
# Toggle this flag: 'solution' or 'error'
PLOT_MODE = 'error' 

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

if PLOT_MODE == 'error':
    # Plot Error Surface
    surf = ax.plot_surface(X, Y, U_error, cmap=cm.inferno, linewidth=0, antialiased=False)
    ax.set_title(f"Error Surface (|BEM - Exact|)\nMax Error: {np.max(U_error):.2e}")
    ax.set_zlabel("Absolute Error")
else:
    # Plot Solution Surface
    surf = ax.plot_surface(X, Y, U_bem, cmap=cm.viridis, linewidth=0, antialiased=False)
    ax.set_title("BEM Solution Surface (u)")
    ax.set_zlabel("Potential u")

ax.set_xlabel("X")
ax.set_ylabel("Y")
fig.colorbar(surf, shrink=0.5, aspect=10)

plt.tight_layout()
plt.show()
