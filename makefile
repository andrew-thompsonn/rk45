

exe_three_body: main.o util.o optimizer.o integrator.o equations.o 
	gcc -Wall -O3 -o exe_three_body main.o util.o optimizer.o integrator.o equations.o -lm
	rm *.o

main.o: src/main.c src/util.h src/optimizer.h src/integrator.h src/equations.h
	gcc -Wall -O3 -c src/main.c

util.o: src/util.c src/integrator.h src/equations.h
	gcc -Wall -O3 -c src/util.c

optimizer.o: src/optimizer.c src/util.h src/integrator.h
	gcc -Wall -O3 -c src/optimizer.c

integrator.o: src/integrator.c 
	gcc -Wall -O3 -c src/integrator.c

equations.o: src/equations.c src/definitions.h
	gcc -Wall -O3 -c src/equations.c

.PHONY: clean
clean:
	rm exe_three_body
