EPS = inverse_condition.eps \
	  filter.eps \
	  pre_fft_distributions.eps \
	  cartoon.eps \
	  cartoon-top.eps

all_eps: $(EPS)

cartoon.eps: cartoon.gpi cartoon.dat
	gnuplot -e "set terminal epscairo size 4,5; set output '$@'" $<

cartoon-top.eps: cartoon-top.gpi cartoon.dat
	gnuplot -e "set terminal epscairo size 4,2.5; set output '$@'" $<

cartoon.dat: cartoon.m
	octave $<

pre_fft_distributions.eps: pre_fft_distributions.py
	python $<

inverse_condition.eps: inverse_condition.py
	python $<

filter.eps: filtercoeffs.py
	python $<

pasa.cls: pasa-template-20170508.zip
	unzip $<

pasa-template-20170508.zip:
	wget -O $@ https://www.cambridge.org/core/services/aop-file-manager/file/58ef9e6e01d33c6306f01bb2

%.png: %.eps
	convert $< $@
