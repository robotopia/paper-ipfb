EPS = inverse_condition.eps \
	  filter.eps \
	  pre_fft_distributions.eps \
	  cartoon.eps

all_eps: $(EPS)

cartoon.eps: cartoon.gpi
	gnuplot -e "set terminal epscairo size 4,5; set output '$@'" $<

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
