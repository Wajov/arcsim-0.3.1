.PHONY: all clean clean-subdirs

all: lib/libtaucs.a

lib/libtaucs.a:
	cd taucs && ./configure module=\!METIS && $(MAKE)
	mkdir -p include/ lib/
	cp taucs/build/*/*.h taucs/src/*.h include/
	cp taucs/lib/*/libtaucs.a lib/

clean: clean-subdirs
	rm -rf include/* lib/*

clean-subdirs:
	rm -f alglib/cpp/src/*.{o,a}
	rm -rf jsoncpp/{bin,buildscons,dist,libs}
	$(MAKE) -C taucs clean
