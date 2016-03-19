def load_board(filename):
    # Fungerer på filer som er formatert som eksempelet nedenfor.
    # feks:
    #006905010970012305020004860503800020000000000080001907054100070207450093060703100
    # Om du vil teste programmet kan du bruke A som gitt under.
    #A = '006905010970012305020004860503800020000000000080001907054100070207450093060703100'
    
    f = open(filename,'r')
    A = f.read()
    B = [[int(A[i+9*j]) for i in range(9)] for j in range(9)]
    return B

def save_board(filename,A):
    B = ''
    for i in range(9):
        for j in range(9):
            B = B + str(A[i][j])
        
    f = open(filename,'w')
    for element in B:
        f.write(element)
    print("Brettet er lagret.")

def print_board(A):
    print("")
    print("    0 1 2   3 4 5   6 7 8  ")
    print("  +-------+-------+-------+")
    for i in range(9):
        print(i,"|", A[i][0], A[i][1], A[i][2],"|", A[i][3], A[i][4], A[i][5], "|", A[i][6], A[i][7], A[i][8], "|")
        if ((i+1) % 3 == 0):
            print("  +-------+-------+-------+")

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

def get_vertical(A,y): # x = [0,...,8]
    B = [0 for i in range(9)]
    for j in range(9):
        B[j] = A[j][y]
    return B
    

def get_input(A,filename):
    print("Skriv 'avbryt' for å avslutte.")
    print("Skriv 'save' for å lagre og avslutte.")
    print("Skriv 'do' for å utføre et trekk.")
    print("Skriv 'undo' for å slette et trekk.")
            
    choise = input()
    if choise == 'save':
        save_board(filename,A)
    if choise == 'save' or choise == 'avbryt':
        return -2,-1,-1
    if choise == 'do' or choise == 'undo':
        x = int(input("Skriv inn x-koordinat(0-8): "))
        y = int(input("Skriv inn y-koordinat(0-8): "))
        nn = -1
    if choise =='do':
        nn = int(input("Skriv inn et tall(1-9): "))
        return nn,x,y
                
    return nn,x,y

def is_finished(A):
    for row in A:
        for element in row:
            if element == 0:
                return 0
    print("\nGratulerer, du vant!\n")
    return 1
def do_move(A,nn,x,y):
    A1 = [[0 for i in range(9)] for j in range(9)]
    for i in range(9):
        for j in range(9):
            A1[i][j] = A[i][j]
    A1[x][y] = nn
    return A1

def undo_move(A,x,y):
    A1 = [[0 for i in range(9)] for j in range(9)]
    for i in range(9):
        for j in range(9):
            A1[i][j] = A[i][j]
    A1[x][y] = 0
    return A1

def main():
    filename = input("Skriv inn navnet på filen med brettet: ")
    A = load_board(filename)
    nn = 0
    while 1:
        print_board(A)
        nn,x,y = get_input(A,filename)
        if nn == -2:
            return
        elif nn == -1:
            A = undo_move(A,x,y)
        else:
            A = do_move(A,nn,x,y)
            if is_finished(A):
                return
 
main()
print('Spillet er avsluttet.')
