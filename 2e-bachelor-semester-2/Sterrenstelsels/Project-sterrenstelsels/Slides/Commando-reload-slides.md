Om de slides te reloaden gebruik:

```bash
cd Project-sterrenstelsels/Slides && \
/Library/TeX/texbin/pdflatex -interaction=nonstopmode my-slides-APOD-temp.tex >/tmp/compile.log 2>&1 && \
cp my-slides-APOD-temp.pdf my-slides-APOD-preview.pdf && \
echo "✓ PDF bijgewerkt"
```
