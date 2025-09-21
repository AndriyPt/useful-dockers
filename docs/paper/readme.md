## Word document generation

```bash
pandoc ~/workspace/project/docs/paper/paper.tex --citeproc --filter pandoc-crossref --from latex --to docx -o ~/paper_new.docx
```