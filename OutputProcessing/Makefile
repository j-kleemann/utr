CPP=g++
CFLAGS=-Wall -Wconversion -Wsign-conversion -O3
ROOTFLAGS=`root-config --cflags --glibs`

all: getHistogram getHistogramZeroDeg histogramToTxt rootToTxt

getHistogram: GetHistogram.cpp
	$(CPP) -o $@ $< $(CFLAGS) $(ROOTFLAGS)
	cp $@ ..
	cp loopGetHistogram.sh ..

getHistogramZeroDeg: GetHistogramZeroDeg.cpp
	$(CPP) -o $@ $< $(CFLAGS) $(ROOTFLAGS)
	cp $@ ..

histogramToTxt: HistogramToTxt.cpp
	$(CPP) -o $@ $< $(CFLAGS) $(ROOTFLAGS)
	cp $@ ..
	cp loopHistogramToTxt.sh ..

rootToTxt: RootToTxt.cpp
	$(CPP) -o $@ $< $(CFLAGS) $(ROOTFLAGS)
	cp $@ ..

.PHONY: all clean

clean:
	rm getHistogram
	rm ../getHistogram
	rm getHistogramZeroDeg
	rm ../getHistogramZeroDeg
	rm ../loopGetHistogram.sh
	rm histogramToTxt
	rm ../histogramToTxt
	rm ../loopHistogramToTxt.sh
	rm rootToTxt
	rm ../rootToTxt
