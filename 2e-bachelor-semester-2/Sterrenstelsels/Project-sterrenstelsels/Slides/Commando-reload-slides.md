Om de slides te reloaden gebruik (voer uit vanuit de repo-root):

```bash
cd Sterrenstelsels/Project-sterrenstelsels/Slides && \
pdflatex -interaction=nonstopmode my-slides-APOD-temp.tex >/tmp/compile.log 2>&1 && \
cp my-slides-APOD-temp.pdf my-slides-APOD-preview.pdf && \
echo "✓ PDF bijgewerkt"
```

**Let op:** Dit veronderstelt dat `pdflatex` in je `$PATH` beschikbaar is (bijv. via MacTeX op macOS, TeX Live op Linux, of MikTeX op Windows).
