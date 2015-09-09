//Makes nonuniform List for x size. Run with > fileName

#include <stdio.h>
main()
{

double start = 500.0;
double end = 600.0;

double listSize = 10000.0;

double steps = (end - start) / listSize;

double output = start;

        for (int i=0; i<listSize; i++) {
                output = output + steps;
                printf("%f", output);
                printf("\n");
        }
}
