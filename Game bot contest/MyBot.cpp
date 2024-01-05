#include "Othello.h"
#include "OthelloBoard.h"
#include "OthelloPlayer.h"
#include <cstdlib>
#include <ctime>
#include <list>
#include <vector>
#include <algorithm>
#include <random>
#include <functional>
#include <bits/stdc++.h>
#include <unordered_map>
#include <string>
using namespace std;
using namespace Desdemona;


#define INF 1e18
#define hashfEXACT 0
#define hashfALPHA 1
#define hashfBETA 2
#define valUNKNOWN 0

class HashValue {
  public:
    int depth;   
    int flags;
    int value;

    HashValue() {
        depth = -1;
        flags = -1;
        value = -1;
    }

    HashValue(int dep,int flag,int val){
        depth = dep;
        flags = flag;
        value = val;
    }
};

Turn myTurn;

clock_t start,finish;
OthelloBoard globalBoard;
int depthLimit = 6;
int moveNum=0;
bool quickExit=false;

// These were used for debugging purposes
int cachedEval = 0;
int totalEval = 0;
int freshEval = 0;

int cachedAlphaBetaCall = 0;
int totalAlphaBetaCall = 0;
int freshAlphaBetaCall = 0;

// int levelNodeVisitCount[6]={0,0,0,0,0,0};
// ***************************************

unordered_map<string, double> boardStrEvalMap;
unordered_map<string, HashValue> transpositionTable;

// This board representation is used as the key for the hash table and also for board evaluation
string createGridRepresentationOfBoard(char grid[][8],OthelloBoard board,Turn turn,char myTurnflag,char opponentTurnFlag){
    string boardRep="";
    // Mark a board grid with 'm' for my turn, 'y' for opponent's turn and 'e' for empty. This grid will be later used by board evaluation function
    for(int i=0;i<8;i++) {
        for(int j=0;j<8;j++) {
        Coin findTurn = board.get(i,j);
        if(findTurn == turn) grid[i][j] = myTurnflag;
        else if(findTurn == other(turn)) grid[i][j] = opponentTurnFlag;
        else grid[i][j] = 'e';
        boardRep += grid[i][j];
        }
    }
    return boardRep;
}

// Return a value from the hash when the same board position is revisited in a future depth iteration with the *same lookahead depth*.
int ProbeHash(OthelloBoard board,Turn turn,int level, int alpha, int beta){
    int depth = depthLimit - level;
    string boardRep = "";
    char grid[8][8];
    if (level&1) {
        boardRep = createGridRepresentationOfBoard(grid,board,turn,'m','y');
    } else {
        boardRep = createGridRepresentationOfBoard(grid,board,turn,'y','m');
    }

    if((transpositionTable.find(boardRep) == transpositionTable.end())){
        return valUNKNOWN;
    }

    HashValue hashedObj = transpositionTable[boardRep];
    // if( moveNum == 1 && cachedAlphaBetaCall < 5){
    //     cout << "Board rep: " << boardRep << " hashed object: " << hashedObj.depth << " , " << hashedObj.value;
    // }
    if (hashedObj.depth >= depth) {
        if (hashedObj.flags == hashfEXACT){
            return hashedObj.value;
        }
        if ((hashedObj.flags == hashfALPHA) && (hashedObj.value <= alpha)){
            return alpha;
        }
        if ((hashedObj.flags == hashfBETA) && (hashedObj.value >= beta)){
            return beta;
        }
    }
    return valUNKNOWN;

}

// Records the hash with key as board representation string and value as the HashValue object
void RecordHash(OthelloBoard board,Turn turn,int level, int val, int hashf){
    char grid[8][8];
    string boardRep = "";
    if (level&1) {
        boardRep = createGridRepresentationOfBoard(grid,board,turn,'m','y');
    } else {
        boardRep = createGridRepresentationOfBoard(grid,board,turn,'y','m');
    }
    HashValue hvalue((depthLimit - level),hashf,val);
    transpositionTable.insert(pair<string, HashValue>(boardRep, hvalue ));
}

bool canMove(char self, char opp, char *str)  {
	if (str[0] != opp) return false;
	for (int ctr = 1; ctr < 8; ctr++) {
		if (str[ctr] == 'e') return false;
		if (str[ctr] == self) return true;
	}
	return false;
}

bool isLegalMove(char self, char opp, char grid[8][8], int startx, int starty)   {
	if (grid[startx][starty] != 'e') return false;
	char str[10];
	int x, y, dx, dy, ctr;
	for (dy = -1; dy <= 1; dy++)
		for (dx = -1; dx <= 1; dx++)    {
	        // keep going if both velocities are zero
			if (!dy && !dx) continue;
			str[0] = '\0';
			for (ctr = 1; ctr < 8; ctr++)   {
				x = startx + ctr*dx;
				y = starty + ctr*dy;
				if (x >= 0 && y >= 0 && x<8 && y<8) str[ctr-1] = grid[x][y];
				else str[ctr-1] = 0;
			}
			if (canMove(self, opp, str)) return true;
		}
	return false;
}

int numValidMoves(char self, char opp, char grid[8][8])   {
	int count = 0, i, j;
	for(i = 0; i < 8; i++) for(j = 0; j < 8; j++) if(isLegalMove(self, opp, grid, i, j)) count++;
	return count;
}

// It considers several factors to evaluate the board position at the horizon(leaf)
double othelloBoardEvaluator(char grid[8][8],string boardRep)  {
    totalEval++;
    if(!(boardStrEvalMap.find(boardRep) == boardStrEvalMap.end())){
        cachedEval++;
        return boardStrEvalMap[boardRep];
    }

	char my_color = 'm',opp_color = 'y';
    int myTiles = 0, oppTiles = 0, i, j, k, myFrontTiles = 0, oppFrontTiles = 0, x, y;
    double parity = 0.0, corner_captured = 0.0, corner_closeness = 0.0, mobility = 0.0, frontier = 0.0, piece_difference = 0.0;

    int X1[] = {-1, -1, 0, 1, 1, 1, 0, -1};
    int Y1[] = {0, 1, 1, 1, 0, -1, -1, -1};
    // Positional strategy matrix
    int V[8][8] = 	{ { 20, -3, 11, 8, 8, 11, -3, 20 },
    				{ -3, -7, -4, 1, 1, -4, -7, -3 },
    				{ 11, -4, 2, 2, 2, 2, -4, 11 },
    				{ 8, 1, 2, -3, -3, 2, 1, 8 },
    				{ 8, 1, 2, -3, -3, 2, 1, 8 },
    				{ 11, -4, 2, 2, 2, 2, -4, 11 },
    				{ -3, -7, -4, 1, 1, -4, -7, -3 },
    				{ 20, -3, 11, 8, 8, 11, -3, 20 } };


	// Piece difference
    for(i = 0; i < 8; i++)
        for(j = 0; j < 8; j++)  {
            if(grid[i][j] == my_color)  {
                piece_difference += V[i][j];
                myTiles++;
            } 
            else if(grid[i][j] == opp_color)  {
                piece_difference -= V[i][j];
                oppTiles++;
            }
            if(grid[i][j] != 'e')   {
                for(k = 0; k < 8; k++)  {
                    x = i + X1[k]; y = j + Y1[k];
                    if(x >= 0 && x < 8 && y >= 0 && y < 8 && grid[x][y] == 'e') {
                        if(grid[i][j] == my_color)  myFrontTiles++;
                        else oppFrontTiles++;
                        break;
                    }
                }
            }
        }

    // Coin parity
    if(myTiles > oppTiles) parity = (100.0 * myTiles)/(myTiles + oppTiles);
    else if(myTiles < oppTiles) parity = -(100.0 * oppTiles)/(myTiles + oppTiles);

    // Frontiers
    if(myFrontTiles > oppFrontTiles) frontier = -(100.0 * myFrontTiles)/(myFrontTiles + oppFrontTiles);
    else if(myFrontTiles < oppFrontTiles) frontier = (100.0 * oppFrontTiles)/(myFrontTiles + oppFrontTiles);

    // Corner occupancy
    myTiles = oppTiles = 0;
    if(grid[0][0] == my_color) myTiles++;
    else if(grid[0][0] == opp_color) oppTiles++;
    if(grid[0][7] == my_color) myTiles++;
    else if(grid[0][7] == opp_color) oppTiles++;
    if(grid[7][0] == my_color) myTiles++;
    else if(grid[7][0] == opp_color) oppTiles++;
    if(grid[7][7] == my_color) myTiles++;
    else if(grid[7][7] == opp_color) oppTiles++;
    corner_captured = 25 * (myTiles - oppTiles);

    // Corner closeness
    myTiles = oppTiles = 0;
    if(grid[0][0] == 'e')   {
        if(grid[0][1] == my_color) myTiles++;
        else if(grid[0][1] == opp_color) oppTiles++;
        if(grid[1][1] == my_color) myTiles++;
        else if(grid[1][1] == opp_color) oppTiles++;
        if(grid[1][0] == my_color) myTiles++;
        else if(grid[1][0] == opp_color) oppTiles++;
    }
    if(grid[0][7] == 'e')   {
        if(grid[0][6] == my_color) myTiles++;
        else if(grid[0][6] == opp_color) oppTiles++;
        if(grid[1][6] == my_color) myTiles++;
        else if(grid[1][6] == opp_color) oppTiles++;
        if(grid[1][7] == my_color) myTiles++;
        else if(grid[1][7] == opp_color) oppTiles++;
    }
    if(grid[7][0] == 'e')   {
        if(grid[7][1] == my_color) myTiles++;
        else if(grid[7][1] == opp_color) oppTiles++;
        if(grid[6][1] == my_color) myTiles++;
        else if(grid[6][1] == opp_color) oppTiles++;
        if(grid[6][0] == my_color) myTiles++;
        else if(grid[6][0] == opp_color) oppTiles++;
    }
    if(grid[7][7] == 'e')   {
        if(grid[6][7] == my_color) myTiles++;
        else if(grid[6][7] == opp_color) oppTiles++;
        if(grid[6][6] == my_color) myTiles++;
        else if(grid[6][6] == opp_color) oppTiles++;
        if(grid[7][6] == my_color) myTiles++;
        else if(grid[7][6] == opp_color) oppTiles++;
    }
    corner_closeness = -10 * (myTiles - oppTiles);

    // Mobility
    myTiles = numValidMoves(my_color, opp_color, grid);
    oppTiles = numValidMoves(opp_color, my_color, grid);
    if(myTiles > oppTiles) mobility = (100.0 * myTiles)/(myTiles + oppTiles);
    else if(myTiles < oppTiles) mobility = -(100.0 * oppTiles)/(myTiles + oppTiles);

    // final weighted score
    double score = (11 * parity) + (850.724 * corner_captured) + (382.026 * corner_closeness) + (86.922 * mobility) + (78.396 * frontier) + (10 * piece_difference);
    freshEval++;
    // boardStrEvalMap[boardRep]=score;
    // if(!(boardStrEvalMap.find(boardRep) == boardStrEvalMap.end()) && score != boardStrEvalMap[boardRep]){
    //     cout << "BBBBBBBBBUUUUUUUUUUUGGGGGGGGG";
    // }

    boardStrEvalMap.insert(pair<string, double>(boardRep, score ));
    return score;
}

// Used to compare board positions when sorting moves
double compareEvaluator(OthelloBoard board,Turn turn) {
    char grid[8][8];
    string boardRep=createGridRepresentationOfBoard(grid,board,turn,'m','y');
    return othelloBoardEvaluator(grid,boardRep);
}

// Used to sort moves based on board positions after the move was made
bool compare(Move a, Move b) {
    OthelloBoard boardOne = globalBoard,boardTwo = globalBoard;
    boardOne.makeMove(myTurn,a);
    boardTwo.makeMove(myTurn,b);
    return compareEvaluator(boardOne,myTurn)>compareEvaluator(boardTwo,myTurn);
}

// The heart of the code. Main algorithm
double alphaBetaMinMaxWithIterativeDeepeningAndTranspositionTable(OthelloBoard board, Move move, Turn turn, short level, double alpha, double beta) {
    // levelNodeVisitCount[level - 1]++;
    totalAlphaBetaCall++;
    int hashf = hashfALPHA;
    int probedValue = ProbeHash(board,turn,level, alpha, beta);
    if (probedValue != valUNKNOWN){
        cachedAlphaBetaCall++;
        return probedValue;
    }
    freshAlphaBetaCall++;
    finish = clock();
    // When timeout occurs quickly exit all alpha-beta calls and return best move found so far in previous depth iterations
    if(((double)(finish-start)/CLOCKS_PER_SEC)>1.90) {
        // cout << "\n =========MYBOT==============Panic exit level:"<<level << "--at time:--" << ((double)(finish-start)/CLOCKS_PER_SEC);
        quickExit=true;
        if(level&1) return -INF;
        return INF;
    }

    // When the depth is reached, return the heuristic score of the current board position
	if(level == depthLimit) {
		char grid[8][8];
        string boardRep = "";
        if (level&1) {
            boardRep = createGridRepresentationOfBoard(grid,board,turn,'m','y');
        } else {
            boardRep = createGridRepresentationOfBoard(grid,board,turn,'y','m');
        }
        double returnValue = othelloBoardEvaluator(grid,boardRep);
        // Record the value in hash with the flag and depth to be reused later
        RecordHash(board,turn,level, returnValue, hashfEXACT);
		return returnValue;
	}
	board.makeMove(turn,move);
	turn = other(turn);
	list<Move> newMoves = board.getValidMoves(turn);
	list<Move>::iterator iter = newMoves.begin();
	double ret = -INF;
	if(level&1) ret *= -1;
	if(!(newMoves.size())) return ret;
	for(;iter!=newMoves.end() && (!quickExit);iter++) {
		double curr = alphaBetaMinMaxWithIterativeDeepeningAndTranspositionTable(board,*iter,turn,level+1,alpha,beta);
		if(level&1) { // MIN Node
			ret = min(ret,curr);
			beta = min(beta,ret);
            RecordHash(board,turn,level, beta, hashfBETA);
		}
		else { // MAX Node
			ret = max(ret,curr);
			alpha = max(alpha,ret);		
            hashf = hashfEXACT;
		}
		if(beta<=alpha) break;
	}
    RecordHash(board,turn,level, alpha, hashf);
    // if (probedValue != valUNKNOWN && probedValue != ret ){
    //     cout << "BUGGGGGGGGGGGGGGGGGGGGGGG";
    // }

	return ret; 
}

class MyBot: public OthelloPlayer
{
    public:
        /**
         * Initialisation routines here
         * This could do anything from open up a cache of "best moves" to
         * spawning a background processing thread. 
         */
        MyBot( Turn turn );

        /**
         * Play something 
         */
        virtual Move play( const OthelloBoard& board );
    private:
};

MyBot::MyBot( Turn turn )
    : OthelloPlayer( turn )
{
}

Move MyBot::play( const OthelloBoard& board )
{
    start = clock();
    quickExit=false;
    // These were used for debugging purposes
    moveNum++;
    totalEval = 0;
    cachedEval=0;
    freshEval = 0;
    totalAlphaBetaCall=0;
    cachedAlphaBetaCall=0;
    freshAlphaBetaCall=0;
    // cout << "\n Movenum: "<< moveNum;
    // ****************************

    list<Move> moves = board.getValidMoves( turn );
    myTurn = turn;
    globalBoard = board;
    moves.sort(compare);
    list<Move>::iterator it = moves.begin();
    Move bestMove((*it).x,(*it).y);
    double retVal = -INF;
    // Increase depth at each iteration so that we examine a whole sub-tree and obtain best move for that depth
    for(depthLimit=6; (!quickExit) && depthLimit < 20;depthLimit=depthLimit + 2){
        // cout << "\n depthLimit:" << depthLimit;
        double MAX = INF, MIN = -INF;
        OthelloBoard copyBoard = board;
        short level = 1;
        for(;it!=moves.end() && (!quickExit);it++) {
            double currValue = alphaBetaMinMaxWithIterativeDeepeningAndTranspositionTable(copyBoard,*it,turn,level,MIN,MAX);
            if(currValue > retVal) {
                retVal = currValue;
                bestMove = *it;
            }
            copyBoard = board;
        }

        finish = clock();
        if(((double)(finish-start)/CLOCKS_PER_SEC)>1.95) {
            // cout << "\n =========MYBOT==============Panic exit level:"<<level << "--at time:--" << ((double)(finish-start)/CLOCKS_PER_SEC);
            // cout << "\n ==========MYBOT========= Total Eval: "<< totalEval;
            // cout << "\n Cached Eval: "<< cachedEval;
            // cout << "\n Fresh Eval: "<< freshEval;
            // cout << "\n ==========MYBOT========= Total ABCall: "<< totalAlphaBetaCall;
            // cout << "\n ==========MYBOT========= Exit time: " << ((double)(finish-start)/CLOCKS_PER_SEC);
            // cout << "\n Cached ABCall: "<< cachedAlphaBetaCall;
            // cout << "\n Fresh ABCall: "<< freshAlphaBetaCall;
            quickExit=true;
            return bestMove;
        }

        it = moves.begin();
    }
    // cout << "\n ==========MYBOT========= Total Eval: "<< totalEval;
    // cout << "\n Cached Eval: "<< cachedEval;
    // cout << "\n Fresh Eval: "<< freshEval;
    // cout << "\n ==========MYBOT========= Total ABCall: "<< totalAlphaBetaCall;
    // cout << "\n ==========MYBOT========= Exit time: " << ((double)(finish-start)/CLOCKS_PER_SEC);
    // cout << "\n Cached ABCall: "<< cachedAlphaBetaCall;
    // cout << "\n Exit time: " << ((double)(finish-start)/CLOCKS_PER_SEC);
    // cout << "\n Fresh ABCall: "<< freshAlphaBetaCall;
    // cout << "\n ==========MYBOT=============== Level count comparision";
    // for( int f=0;f<6;f++){
    //     cout << "\n Level: " << f << " : " << levelNodeVisitCount[f];
    // }

    return bestMove;
}

// The following lines are _very_ important to create a bot module for Desdemona

extern "C" {
    OthelloPlayer* createBot( Turn turn )
    {
        return new MyBot( turn );
    }

    void destroyBot( OthelloPlayer* bot )
    {
        delete bot;
    }
}