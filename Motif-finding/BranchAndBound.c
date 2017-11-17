#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <stdlib.h>

#define k 4 //alphabet size
#define INFINITE 10000 //reasonably big number. No significance.

void nextLeaf(int *a, int L)
{
	for ( int i = (L-1); i >= 0; --i)
	{
		if (a[i] < k)
		{
			a[i] = a[i] + 1;
			return;
		}

		a[i] = 1;
	}
}

void allLeaves(int L)
{
	int a[L];
	int allDone = 0;
	for (int i = 0; i < L; ++i)
	{
		a[i] = 1;
	}

	while(1)
	{
		allDone = 1;
		nextLeaf(a, L);

		printf( "Vertex: ");
		for( int i = 0; i < L; ++i)
		{
			printf( "%d ", a[i]);
			if (a[i] != 1)
			{
				allDone = 0;
			}
		}
		printf( "\n");

		if( allDone == 1)
			return;
	}
}

void printVertex(int *a, int L)
{
	printf("Vertex: ");

	for ( int i = 0; i < L; ++i)
	{
		printf("%d ", a[i]);
	}

	printf( "\n");
}

void nextVertex( int *a, int *i, int L)
{
	if (*i < L)
	{
		a[*i] = 1;
		*i = *i + 1;
		return;
	}
	else
	{
		for (int j = (L -1); j >= 0; --j)
		{
			if (a[j] < k)
			{
				a[j] = a[j] + 1;
				*i = j+1;
				return;
			}
		}
	}
	return;
}

void bypass(int *a, int *i, int L)
{
	for (int j = (*i-1); j >= 0; --j)
	{
		if (a[j] < k)
		{
			a[j] = a[j] + 1;
			*i = j+1;
			return;
		}
	}

	*i = 0;
}

int hammingDistance(int* v, int *t, int L, int seq_size)
{
	int score = 0;
	int bestScore = 0;

	for( int i = 0; i < (seq_size - L); ++i)
	{
		for( int j = 0; j < L; ++j)
		{
			if ( v[j] != t[i+j])
			{
				++score;
			}
		}

		if ( i == 0)
		{
			bestScore = score;
		}
		else if ( score < bestScore)
		{
			bestScore = score;
		}

		score = 0;
	}

	return bestScore;
}

int totalHammingDistance( int *v, int** dna, int dna_size, int L, int seq_size)
{
	int score = 0;

	for ( int i = 0; i < dna_size; ++i)
	{
		score = score + hammingDistance(v, dna[i], L, seq_size);
	}

	return score;
}

void buildBestArray (int* bestArray, int *lmer, int *level, int L, int half, int **dna, int dna_size, int seq_size)
{
	int currentScore = INFINITE;
	
	if ( *level == 0)
	{
		return;
	}
	else
	{
		if (  (*level+1) == half)
		{
			currentScore = totalHammingDistance(lmer, dna, dna_size, *level, seq_size);
			if ( bestArray[*level] > currentScore)
			{
				bestArray[*level] = currentScore;
			}
			
			if ( lmer[*level-1] == k)
			{
				*level = *level - 1;
				bypass(lmer, level, L);
			}
			else
			{
				bypass(lmer, level, L);
			}
			buildBestArray(bestArray, lmer, level, L, half, dna, dna_size, seq_size);
		}
		else if ( *level < half)
		{
			currentScore = totalHammingDistance(lmer, dna, dna_size, *level, seq_size);
			if ( bestArray[*level] > currentScore)
			{
				bestArray[*level] = currentScore;
			}

			nextVertex(lmer, level, L);
			buildBestArray(bestArray, lmer, level, L, half, dna, dna_size, seq_size);
		}
		if ( *level >= half)
		{
			*level = 0;
		}
		
		return;
	}
}

int medianStringSearch(int **dna, int* bestWord, int L, int dna_size, int seq_size)
{
	int lmer[L];
	int bestSubstring[L];
	int bestDistance = INFINITE;
	int level = 1;

	int half;

	if ( (L%2) != 0)
	{
		half = (L/2) + 1;
	}
	else
	{
		half = (L/2);
	}

	for (int i = 0; i < L; ++i)
	{
		lmer[i] = 1;
		bestWord[i] = 0;
		bestSubstring[i] = INFINITE;
	}

	buildBestArray(bestSubstring, lmer, &level, L, half, dna, dna_size, seq_size);

	for ( int i = 0; i < L; ++i)
	{
		lmer[i] = 1;
	}

	int optimisticPrefixDistance, optimisticSuffixDistance, newDistance;
	level = half;
	while (level > 0)
	{
		if ( (level < L) )
		{
			int prefix[level];

			//prefix
			for ( int i = 0; i < level; ++i)
			{
				prefix[i] = lmer[i];
			}
			
			optimisticPrefixDistance = totalHammingDistance(lmer, dna, dna_size, level, seq_size);	

			if ( optimisticPrefixDistance < bestSubstring[level] )
			{
				bestSubstring[level] = optimisticPrefixDistance;
			}

			if ( (L-level) < level)
			{
				optimisticSuffixDistance = bestSubstring[L-level];
			}
			else
			{
				optimisticSuffixDistance = 0;
			}

			if ( (optimisticPrefixDistance + optimisticSuffixDistance) >= bestDistance)
			{
				bypass(lmer, &level, L);
			}
			else
			{
				nextVertex(lmer, &level, L);
			}
		}
		else
		{
			newDistance = totalHammingDistance( lmer, dna, dna_size, L, seq_size);

			if ( bestDistance > newDistance)
			{
				bestDistance = newDistance;

				for ( int i = 0; i < L; ++i)
				{
					bestWord[i] = lmer[i];
				}
			}

			nextVertex(lmer, &level, L);
		}
	}

	return bestDistance;

}

int lettersToNumbers( char c)
{
	if ( c == 'A') return 1;
	if ( c == 'C') return 2;
	if ( c == 'G') return 3;
	if ( c == 'T') return 4;
}

char numberToLetter( int c)
{
	if ( c == 1 ) return 'A';
	if ( c == 2 ) return 'C';
	if ( c == 3 ) return 'G';
	if ( c == 4 ) return 'T';
}
