#!/usr/bin/env python3

def needleman_wunsch(seq_a: str, seq_b: str, match: int, mismatch: int, gap: int) -> tuple[tuple[str, str], int]:
    matrix = generate_matrix(len(seq_a), len(seq_b))
    filled_matrix = fill_matrix(matrix = matrix, seq_a = seq_a, seq_b = seq_b, match = match, mismatch = mismatch, gap = gap)
    aln, score = find_optimal_alignment(matrix = filled_matrix, seq_a = seq_a, seq_b = seq_b, match = match, mismatch = mismatch, gap = gap)

    # ret_tup = (('CTTCTCGT-CGGTCTCGTGGTTCGGGAAC', 'CTT-TCATCCACT-TCGTTGCCCGGGAAC'), 11)
    return aln, score

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
        matrix[0][col] = matrix[0][col - 1] + gap

    for row in range(1, len(seq_a) + 1):
        matrix[row][0] = matrix[row - 1][0] + gap

    for row in range(1, len(seq_a) + 1):
        seq_a_base = seq_a[row - 1]
        for col in range(1, len(seq_b) + 1):
            seq_b_base = seq_b[col - 1]
            matrix[row][col] = max(
                                    (matrix[row - 1][col] + gap), 
                                    (matrix[row][col - 1] + gap), 
                                    (matrix[row - 1][col -1] + match if seq_a_base == seq_b_base else matrix[row - 1][col -1] + mismatch )
                                    )    
    # for row in matrix:
    #     print(row)    
    return matrix

def find_optimal_alignment(matrix: list[list[list[int]]], seq_a: str, seq_b: str, match: int, mismatch: int, gap: int) -> tuple[tuple[str, str], int]:
    score = matrix[len(seq_a)][len(seq_b)]
    # print(score)

    aln_a = ""
    aln_b = ""
    row = len(seq_a)
    col = len(seq_b)
    while row > 0:
        while col > 0:
            l = [matrix[row-1][col-1], matrix[row][col-1], matrix[row-1][col]]
            # max_val = max(l)
            max_pos = l.index(max(l))
            #match
            if max_pos == 0 and matrix[row][col] == matrix[row-1][col-1] + match:
                aln_b += seq_b[col-1]
                aln_a += seq_a[row-1]
                row -= 1
                col -= 1
            #mismatch
            elif matrix[row][col] == matrix[row-1][col-1] + mismatch:
                aln_b += seq_b[col-1]
                aln_a += seq_a[row-1]
                row -= 1
                col -= 1
            #gap on seq_b
            elif max_pos == 1 and matrix[row][col] == matrix[row][col-1] + gap:
                aln_a += "-"
                aln_b += seq_b[col-1]
                col -= 1
            #gap on seq_a
            elif max_pos == 2 and matrix[row][col] == matrix[row-1][col] + gap:
                aln_b += "-"
                aln_a += seq_a[row-1]
                row -= 1
            #deafult to moving diagonally up
            else:
                aln_b += seq_b[col-1]
                aln_a += seq_a[row-1]
                row -= 1
                col -= 1
    # print(f"{aln_a}\n{aln_b}")
    a = aln_a[::-1]
    b = aln_b[::-1]
    return ((a, b), score)

# if __name__ == "__main__":
#     needleman_wunsch("ATGA", "ATCA", 1, -1, -1)
#     needleman_wunsch("TAGTCAT", "TATCAAT", 1, -1, -1)
#     needleman_wunsch("CTTCTCGTCGGTCTCGTGGTTCGGGAAC", "CTTTCATCCACTTCGTTGCCCGGGAAC", 1, -1, -1)
