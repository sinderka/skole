def load_board(filename):

    # Works on files formated as shown below.
    # Ex:
    # 600020040300007100051000706000005080200000001040300000406000970002600008090070004
    # The program can be tested by using, if files are not desired
    # A = '600020040300007100051000706000005080200000001040300000406000970002600008090070004'
    f = open(filename, 'r')
    A = f.read()
    B = [[int(A[i+9*j]) for i in range(9)] for j in range(9)]
    return B


def print_board(A):
    print("")
    print("    0 1 2   3 4 5   6 7 8  ")
    print("  +-------+-------+-------+")
    for i in range(9):
        print(i,"|", A[i][0], A[i][1], A[i][2],"|", A[i][3], A[i][4], A[i][5], "|", A[i][6], A[i][7], A[i][8], "|")
        if ((i+1) % 3 == 0):
            print("  +-------+-------+-------+")


def print_moves(nn,x,y):
    print("nn: ",nn,"x: ",x,"y: ", y)


def get_square(A,x,y): # x = [0,1,2], y = [0,1,2]
    B = [0 for j in range(9)]
    for i in range(3):
        for j in range(3):
                B[i+3*j] = A[x*3+i][y*3+j]
    return B


def get_horizontal(A,x): # x = [0,...,8]
    B = [0 for i in range(9)]
    for i in range(9):
        B[i] = A[x][i]
    return B


def get_vertical(A,y): # y = [0,...,8]
    B = [0 for i in range(9)]
    for j in range(9):
        B[j] = A[j][y]
    return B


def is_legal(A,nn,x,y):
    if (A[x][y] != 0) or (nn in get_horizontal(A,x)) or (nn in get_vertical(A,y)) or (nn in get_square(A,round(x/3-0.4),round(y/3-0.4))):
        return 0
    return 1


def is_finished(A):
    for row in A:
        for element in row:
            if element == 0:
                return 0
    print_board(A)
    print("\nGratulerer, du vant!\n")
    return 1

def mat_copy(D,A,dept):
    D1 = [[[0 for k in range(9)] for i in range(9)] for j in range(dept+1)]
    for k in range(dept):
        for i in range(9):
            for j in range(9):
                D1[k][i][j] = D[k][i][j]
    for i in range(9):
        for j in range(9):
            D1[dept][i][j] = A[i][j]
    return D1

def do_move(A,nn,x,y):
    A1 = [[0 for i in range(9)] for j in range(9)]
    for i in range(9):
        for j in range(9):
            A1[i][j] = A[i][j]
    A1[x][y] = nn
    return A1

def verify_board(A):
    for i in range(9):
        for j in range(9):
            nn = A[i][j]
            A[i][j] = 0
            if not (nn == 0) and not is_legal(A,nn,i,j):
                print_board(A)
                return 0
            nn, A[i][j] = 0, nn
    return 1

def all_possible_moves(A):
    possible_moves = [[[(k+1)*is_legal(A,k+1,i,j) for k in range(9)] for i in range(9)] for j in range(9)]
    return possible_moves

def fewest_possible_moves(possible_moves,wrong_moves):
    fewest_moves = 1000
    nn,x,y = 0,0,0

    for i in range(9):
        for j in range(9):
            temp = 0
            for k in range(9):

                if possible_moves[i][j][k]>0:
                    temp += 1
            if fewest_moves > temp and temp >= 1:
                fewest_moves = temp
                nn,x,y = max(possible_moves[i][j][:]),j,i
                while wrong_moves != -1 and ([nn,x,y] in wrong_moves) :
                    possible_moves[i][j][possible_moves[i][j][:].index(nn)] = 0
                    nn,x,y = max(possible_moves[i][j][:]),j,i

    return nn,x,y,fewest_moves

def update_last_moves(last_moves,nn,x,y,dept):
    if dept == 1:
        return [[nn,x,y]]
    last_moves1 = [[0 for i in range(3)] for j in range(dept)]
    for i in range(dept-1):
        last_moves1[i][0],last_moves1[i][1],last_moves1[i][2] = last_moves[i][0],last_moves[i][1],last_moves[i][2]
    last_moves1[dept-1][0],last_moves1[dept-1][1],last_moves1[dept-1][2]= nn,x,y
    return last_moves1


def next_move(A,possible_moves,dept,last_moves,wrong_moves):
    try:
        nn,x,y,fewest_moves = fewest_possible_moves(possible_moves,wrong_moves[dept])
    except:
        nn,x,y,fewest_moves = fewest_possible_moves(possible_moves,-1)

    if nn == 0:
        return last_moves,A,dept,1
    A = do_move(A,nn,x,y)

    dept += 1

    last_moves =update_last_moves(last_moves,nn,x,y,dept)

    if verify_board(A):
        return last_moves,A,dept,0
    print("arg!!!!")

def main(filename):
    A = load_board(filename)
    print_board(A)

    dept = 0


    wrong_moves = {}
    last_moves=[[]]

    while is_finished(A) == 0:

        possible_moves = all_possible_moves(A)
        last_moves,A,dept,back = next_move(A,possible_moves,dept,last_moves,wrong_moves)
        print('dept',dept)
        if back:
            print('backtracking')
            nn,x,y = last_moves[:][dept-1]
            A = do_move(A,0,x,y)
            last_moves = last_moves[:][0:dept-1]
            dept -=1


            try:
                wrong_moves[dept].append([nn,x,y])
            except:
                wrong_moves[dept] = [[nn,x,y]]


            for key in list(wrong_moves.keys()):
                if key > dept:
                    del wrong_moves[key]


main('sudoku3.txt')