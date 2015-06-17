#include <stdio.h>
main()
{

double start = 64.325;
double end = 298.0;

double listSize = 10000.0;

double steps = (end - start) / listSize;

double output = start;

	for (int i=0; i<listSize; i++) {
		output = output + steps;
		printf("%f", output);
		printf("\n");
	}
}
