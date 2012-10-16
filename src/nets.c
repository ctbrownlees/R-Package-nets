
#include <R.h>

#include <stdio.h>
#include <stdlib.h>

void adding(double* a, int* b, double* ab)
{
	int i;
	for (i = 0; i < *b; i++)
	{
		ab[i] = a[i] + 5;
	}
}

