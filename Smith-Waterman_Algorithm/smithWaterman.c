#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#define MIN_SCORE  0
#define MAX(X,Y) ( ( (X) > (Y) ) ? (X) : (Y) )

enum pointer { left, up, diagnol, end};

int subMatrix[4][4] = { { 4, -2, -3, -1}, {-2,  4, -1, -3}, {-3, -1,  4, -2}, {-1, -3, -2,  4} };


struct cell
{
	int score;
	enum pointer backtrackPointer;
};

// Because the alignment can be of empty sub-sequences.
void baseCondition( struct cell **matrix, int rows, int cols)
{
	//initialize first column to zeroes
	int firstCol = 0;
	for ( int i = 0; i < rows; ++i)
	{
		(matrix[i][firstCol]).score = MIN_SCORE;
		(matrix[i][firstCol]).backtrackPointer = end;
	}

	//initialize first row to zeroes;
	int firstRow = 0;
	for ( int i = 0; i < cols; ++i)
	{
		(matrix[firstRow][i]).score = MIN_SCORE;
		(matrix[firstRow][i]).backtrackPointer = end;
	}
}

// Mapping A-0, T-1, C-2, G-3
int lookUpScore( int indexRow, int indexCol)
{
	return subMatrix[ ( (indexRow == 71) ? (3) : ( (indexRow - 65) % 3)) ][ (indexCol == 71) ? (3) : ( (indexCol - 65) % 3) ];
}

// Figure out relevant gap penalty for insertion.
int affineLeft( struct cell* c, int gapStart)
{
	int start = c->score + gapStart;
	int extend = c->score;

	if ( c->backtrackPointer == left)
	{
		return extend;
	}
	else {
		return start;
	}
}

// Figure out relevant gap penalty for deletion.
int affineUp( struct cell* c, int gapStart)
{
	int start = c->score + gapStart;
	int extend = c->score;

	if ( c->backtrackPointer == up) { 
		return extend;
	}
	else
	{
		return start;
	}
}


// V(i, j) = max[ 0, E(i, j), F(i, j), G(i, j)]
void recurV(struct cell **matrix, int row , int col, char naiveOrAffine, int* rowLetter, int* colLetter)
{
	int gapStart = -16;
	int gapExtend = -4;

	int leftScore, upScore, diagScore, subScore;

	//G(i, j) = V(i-1, j-1) + score(i-1, j-1)
	subScore = lookUpScore( rowLetter[row], colLetter[col] ); //look up match/mismatch score in the substitution matrix.
	diagScore = subScore + matrix[row-1][col-1].score;
	
	//////// gap score calculation with affine gap penalty //////////
	if ( naiveOrAffine == 'a') { 
		leftScore = affineLeft( &matrix[row][col-1], gapStart) + gapExtend;
		upScore = affineUp( &matrix[row-1][col], gapStart) + gapExtend;
	}
	////// gap score calculation with naive gap penalty ////////
	else {

		leftScore = matrix[row][col-1].score + gapExtend;

		upScore = matrix[row-1][col].score + gapExtend;
	}
	
	//max[0, diagScore, upScore, leftScore]
	if ( leftScore > MAX(diagScore, upScore) )
	{
		matrix[row][col].score = leftScore;
		matrix[row][col].backtrackPointer = left;
	}
	else if ( upScore > MAX(diagScore, leftScore))
	{
		matrix[row][col].score = upScore;
		matrix[row][col].backtrackPointer = up;
	}
	else
	{
			matrix[row][col].score = diagScore;
			matrix[row][col].backtrackPointer = diagnol;
			if ( (upScore == leftScore) && upScore > diagScore)
			{
				matrix[row][col].score = upScore;
				matrix[row][col].backtrackPointer = up;
			}
	}

	if ( matrix[row][col].score <= MIN_SCORE)
	{
		matrix[row][col].score = MIN_SCORE;
		matrix[row][col].backtrackPointer = end;
	}
	
}

//BackTrack recursively.
//A pointer to the left corresponds to a dash in the first string.
//A pointer upwards corresponds to a dash in the second string.
//A diagnol pointer means either a match or a mismatch.
void traceBack( struct cell **matrix, int i, int j, int* rowLetter, int* colLetter, int* colTrack, int* trackLength)
{
	if ( matrix[i][j].backtrackPointer == end)
	{
		return;
	}
	else if ( matrix[i][j].backtrackPointer == diagnol)
	{
		colTrack[++(*trackLength)] = colLetter[j];
		traceBack( matrix, i-1, j-1, rowLetter, colLetter, colTrack, trackLength);
		printf( "%c", (char) rowLetter[i]);
	}
	else if ( matrix[i][j].backtrackPointer == up)
	{
		colTrack[++(*trackLength)] = '-';
		traceBack( matrix, i-1, j, rowLetter, colLetter, colTrack, trackLength);
		printf( "%c", (char) rowLetter[i] );
	}
	else if ( matrix[i][j].backtrackPointer == left)
	{
		colTrack[++(*trackLength)] = colLetter[j];
		traceBack( matrix, i, j-1, rowLetter, colLetter, colTrack, trackLength);
		printf( "-");
	}
	//debugging
	else
	{
		printf("WTF");
		return;
	}
	
}

// Not really recursive.....
void dynamicRecurrence( struct cell **matrix, int rows, int cols, char naiveOrAffine, int* rowLetter, int* colLetter)
{
	int seqStart = 1;

	//set base condition, i.e V(i,0) = V(0, j) = 0
	baseCondition( matrix, rows, cols);

	//begin filling the matrix row-wise with the following recurrences
	// V(i, j) = max[ E(i, j), F(i, j), G(i, j)]
	// E[i, j) = max[ E(i, j-1)+Ws, V(i, j-1) + Wg] 
	// F(i, j) = max[ F(i-1, j)+Ws, V(i-1, j) + Wg] 
	// G(i, j) = V(i-1, j-1) + (S1(i) == S2(j) ? (W-match) : (W-mismatch)
	
	int r = 0;
	int c = 0;
	for ( int i = seqStart; i < rows; ++i)
	{
		for ( int j = seqStart; j < cols; ++j)
		{
			recurV( matrix, i, j, naiveOrAffine, rowLetter, colLetter);
			if ( matrix[r][c].score < matrix[i][j].score)
			{
				r = i;
				c = j;
			}
		}
	}
	//save values of sequence two, i.e the columns, to print later.
	int colTrack[cols];
	int trackLength = -1;

	//print max score.
	printf ( "%dx", matrix[r][c].score);

	traceBack( matrix, r , c, rowLetter, colLetter, colTrack, &trackLength);
	printf("x");
	for ( int i = trackLength; i >= 0; --i)
	{
		printf( "%c", (char) colTrack[i]);
	}
	printf( "x");
}

int main()
{
	//put first sequence in intermediary file seq.txt.
	system( "cat input.fa | tr -d '\n' | grep -o -E \"1[ATCG]+\" | cut -d '1' -f2 > seq.txt");

	//put second sequence in intermediary file seq.txt.
    system( "cat input.fa | tr -d '\n' | grep -o -E \"2[ATCG]+\" | cut -d '2' -f2 >> seq.txt");

    int seq_size1 = 0;
    int seq_size2 = 0;
    FILE *stream;
    struct stat sb;
    size_t fileSize;
    char *buffer;

	///////////////////////////////////////////////
	//											 //
	// process file seq.txt to get relevant data //
	//											 //
	///////////////////////////////////////////////
    if ( stat("seq.txt", &sb) == 0)
    {
            fileSize = sb.st_size;
    }
    else
    {
        printf("Exiting. 'seq.txt' not found");
        exit(1);
    }

    buffer = (char *)malloc(fileSize * sizeof(char));
    if ( buffer == NULL)
    {
        perror( "Unable to allocate buffer.");
        exit(1);
    }

    //process first sequence.
    stream = fopen( "seq.txt", "r");
    seq_size1 = getline(&buffer, &fileSize, stream);

    int row[seq_size1];
    row[0] = -1;

    for ( int i = 1; i < seq_size1; ++i)
    {
        row[i] = buffer[i-1];
    }

    //process second sequence.
    seq_size2 = getline(&buffer, &fileSize, stream);

    fclose(stream);

    int col[seq_size2];
    col[0] = -1;

    for ( int i = 1; i < seq_size2; ++i)
    {
        col[i] = buffer[i-1];
    }

    free(buffer);


	////////////////////////////////
	//							  //
	//	Begin algorithm execution //
	//							  //
	////////////////////////////////
	

	//allocatate mem for matrix.
	/* NOTE: THIS MAY NOT WORK WITH VERY LARGE VALUES OF SEQ_SIZE1 AND SEQ_SIZE2.
	 * -> Malloc doesnt return contiguous memory blocks, therefore for large sizes
	 *  of seq_size1 and seq_size2, we might end up with discontiguous memory segments. The subsequent code
	 *  assumes a continuous memory,and would thus cause segmentation faults.
	 * FIX: USE ONE DIMENSIONAL ARRAY AND RELEVANT ARITHMETIC TO GIVE MATRIX LIKE STRUCTURE.
	 * -> Unfortunately, I am too lazy to implement this :( .
	*/

    struct cell **matrix;
	matrix = malloc( seq_size1 * sizeof(*matrix));

    for ( int i = 0; i < seq_size1; ++i)
    {
        matrix[i] = malloc ( seq_size2 * sizeof(matrix[i]) );
    }

	//The real meat starts here. 
    printf( "Output1");
    dynamicRecurrence( matrix, seq_size1, seq_size2, 'f', row, col);
    printf( "Output2");
    dynamicRecurrence( matrix, seq_size1, seq_size2, 'a', row, col);

	//deallocs, cause CS-202 left us traumatized.
	for ( int i = 0; i < seq_size1; ++i)
	{
		free(matrix[i]);
	}
	free(matrix);
    return 0;
}
