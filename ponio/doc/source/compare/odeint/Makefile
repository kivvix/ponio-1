all: lorenz.txt  transport.txt arenstorf.txt

%.txt: %.out
	./$<

%.out: %.cpp
	${CXX} -std=c++20 $< -I ${CONDA_PREFIX}/include -o $@

clean:
	rm -f *.txt *.out
