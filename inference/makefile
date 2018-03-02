CC = clang++
CFLAGS = -std=c++11 -O3
EXTRA =

all:
	make -j16 all2

s:
	make CC=g++ EXTRA="-pthread" -j16 all2

copy:
	scp makefile *.cpp *.h ajbean@csp-serv-01.csp.illinois.edu:~/inheritance/

all2: GraphSolver.o SumProductBP.o BetaStepBP.o ProjectionBP.o AProjectionBP.o ZoomBP.o FastZoomBP.o KProjectionBP.o StochasticBP.o StochasticBP2.o Energy.bin StepDist.bin Disparities.bin test.bin




DiscreteFactorGraph.o: DiscreteFactorGraph.cpp DiscreteFactorGraph.h FloatType.h
	$(CC) $(CFLAGS) $(EXTRA) -c DiscreteFactorGraph.cpp -o DiscreteFactorGraph.o

GraphSolver.o: GraphSolver.cpp GraphSolver.h DiscreteFactorGraph.h FloatType.h
	$(CC) $(CFLAGS) $(EXTRA) -c GraphSolver.cpp -o GraphSolver.o




SumProductBP.o: SumProductBP.cpp SumProductBP.h DiscreteFactorGraph.h FloatType.h GraphSolver.h
	$(CC) $(CFLAGS) $(EXTRA) -c SumProductBP.cpp -o SumProductBP.o

BetaStepBP.o: BetaStepBP.cpp BetaStepBP.h DiscreteFactorGraph.h FloatType.h GraphSolver.h
	$(CC) $(CFLAGS) $(EXTRA) -c BetaStepBP.cpp -o BetaStepBP.o

ProjectionBP.o: ProjectionBP.cpp ProjectionBP.h DiscreteFactorGraph.h FloatType.h GraphSolver.h
	$(CC) $(CFLAGS) $(EXTRA) -c ProjectionBP.cpp -o ProjectionBP.o

AProjectionBP.o: AProjectionBP.cpp AProjectionBP.h DiscreteFactorGraph.h FloatType.h GraphSolver.h
	$(CC) $(CFLAGS) $(EXTRA) -c AProjectionBP.cpp -o AProjectionBP.o

ZoomBP.o: ZoomBP.cpp ZoomBP.h DiscreteFactorGraph.h FloatType.h GraphSolver.h
	$(CC) $(CFLAGS) $(EXTRA) -c ZoomBP.cpp -o ZoomBP.o

FastZoomBP.o: FastZoomBP.cpp FastZoomBP.h DiscreteFactorGraph.h FloatType.h GraphSolver.h
	$(CC) $(CFLAGS) $(EXTRA) -c FastZoomBP.cpp -o FastZoomBP.o

# NewZoomBP.o: NewZoomBP.cpp NewZoomBP.h DiscreteFactorGraph.h FloatType.h GraphSolver.h
#   $(CC) $(CFLAGS) $(EXTRA) -c NewZoomBP.cpp -o NewZoomBP.o

KProjectionBP.o: KProjectionBP.cpp KProjectionBP.h DiscreteFactorGraph.h FloatType.h GraphSolver.h
	$(CC) $(CFLAGS) $(EXTRA) -c KProjectionBP.cpp -o KProjectionBP.o

StochasticBP.o: StochasticBP.cpp StochasticBP.h DiscreteFactorGraph.h FloatType.h GraphSolver.h
	$(CC) $(CFLAGS) $(EXTRA) -c StochasticBP.cpp -o StochasticBP.o

StochasticBP2.o: StochasticBP2.cpp StochasticBP2.h DiscreteFactorGraph.h FloatType.h GraphSolver.h
	$(CC) $(CFLAGS) $(EXTRA) -c StochasticBP2.cpp -o StochasticBP2.o




Energy.bin: Energy.cpp Solvers.h DiscreteFactorGraph.o GraphSolver.o SumProductBP.o BetaStepBP.o ProjectionBP.o AProjectionBP.o ZoomBP.o FastZoomBP.o KProjectionBP.o StochasticBP.o StochasticBP2.o
	$(CC) $(CFLAGS) $(EXTRA) -o Energy.bin Energy.cpp DiscreteFactorGraph.o GraphSolver.o SumProductBP.o BetaStepBP.o ProjectionBP.o AProjectionBP.o ZoomBP.o FastZoomBP.o KProjectionBP.o StochasticBP.o StochasticBP2.o

StepDist.bin: StepDist.cpp Solvers.h DiscreteFactorGraph.o GraphSolver.o SumProductBP.o
	$(CC) $(CFLAGS) $(EXTRA) -o StepDist.bin StepDist.cpp DiscreteFactorGraph.o GraphSolver.o SumProductBP.o

Disparities.bin: Disparities.cpp Solvers.h DiscreteFactorGraph.o GraphSolver.o SumProductBP.o BetaStepBP.o ProjectionBP.o AProjectionBP.o ZoomBP.o FastZoomBP.o KProjectionBP.o StochasticBP.o StochasticBP2.o
	$(CC) $(CFLAGS) $(EXTRA) -o Disparities.bin Disparities.cpp DiscreteFactorGraph.o GraphSolver.o SumProductBP.o BetaStepBP.o ProjectionBP.o AProjectionBP.o ZoomBP.o FastZoomBP.o KProjectionBP.o StochasticBP.o StochasticBP2.o

test.bin: test.cpp
	$(CC) $(CFLAGS) -o test.bin test.cpp

clean:
	rm -rf *.o *.bin *.dSYM

