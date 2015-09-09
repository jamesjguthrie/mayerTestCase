#include <stdio.h>
main()
{

double start = 400;
double end = 500;

double listSize = 10000.0;

double steps = (end - start) / listSize;

double output = start;

	for (int i=0; i<listSize; i++) {
		output = output + steps;
		printf("%f", output);
		printf("\n");
	}
}
