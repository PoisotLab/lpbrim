PKG=lpbrim

ALL: $(PKG).tar.gz

doc: R/*.r
	R -e "library(devtools); document('.')"

cran/$(PKG):
	mkdir -p cran/$(PKG)

test: cran/$(PKG) doc
	cp -r * cran/$(PKG) 2>/dev/null; true
	rm -r cran/$(PKG)/{cran,tests} 2>/dev/null; true
	rm cran/$(PKG)/Makefile 2>/dev/null; true
	rm cran/$(PKG)/$(PKG).tar.gz 2>/dev/null; true
	cd cran; R CMD check --as-cran $(PKG)

$(PKG).tar.gz:
	rm $@ 2>/dev/null; true
	cd cran; tar -zcvf $@ $(PKG)
	mv cran/$(PKG).tar.gz $@

clean:
	rm $(PKG).tar.gz
	rm -r cran/$(PKG)
