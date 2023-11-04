#!/usr/bin/env python3

def needleman_wunsch(seq_a: str, seq_b: str, match: int, mismatch: int, gap: int) -> tuple[tuple[str, str], int]:
    matrix = generate_matrix(len(seq_a), len(seq_b))
    filled_matrix = fill_matrix(matrix = matrix, seq_a = seq_a, seq_b = seq_b, match = match, mismatch = mismatch, gap = gap)


    # ret_tup = (('CTTCTCGT-CGGTCTCGTGGTTCGGGAAC', 'CTT-TCATCCACT-TCGTTGCCCGGGAAC'), 11)
    # return ret_tup

def generate_matrix(row: int, col: int) -> list[list[list[None]]]:
    """
    creates a (row + 1) X (col + 1) matrix with each cell as a list.
    extra row and column added to fill in initial values.

    Args:
        row: length of seq_a
        col: length of seq_b

    Returns:
        3-level nested matrix
    """
    mat = [[ 0 for j in range(col + 1)] for i in range(row + 1)]
    # print(mat)
    return mat

def fill_matrix(matrix: list[list[list[None]]], seq_a: str, seq_b: str, match: int, mismatch: int, gap: int) -> list[list[list[int]]]:
    for col in range(1, len(seq_b) + 1):
        matrix[0][col] = matrix[0][col - 1] - 1

    for row in range(1, len(seq_a) + 1):
        matrix[row][0] = matrix[row - 1][0] - 1

    for row in range(1, len(seq_a) + 1):
        seq_a_base = seq_a[row - 1]
        for col in range(1, len(seq_b) + 1):
            seq_b_base = seq_b[col - 1]
            matrix[row][col] = max(
                                    (matrix[row - 1][col] - 1), 
                                    (matrix[row][col - 1] - 1), 
                                    (matrix[row - 1][col -1] + 1 if seq_a_base == seq_b_base else matrix[row - 1][col -1] - 1 )
                                    )    
    for row in matrix:
        print(row)    
    return

if __name__ == "__main__":
    needleman_wunsch("ATA", "ATCA", 1, -1, -1)
    print('\n')
    needleman_wunsch("TAGTCAT", "TATCAAT", 1, -1, -1)