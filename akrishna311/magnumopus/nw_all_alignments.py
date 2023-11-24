#!/usr/bin/env python3

def needleman_wunsch(seq_a: str, seq_b: str, match: int, mismatch: int, gap: int) -> tuple[tuple[str, str], int]:
    """
    function to implement NW algorithm

    Args:
        seq_a: first sequence to align
        seq_b: second sequence to align
        match: match score
        mismatch: mismatch penalty
        gap: gap penalty

    Returns:
        A tuple containing aligned sequences and the matched score
    """
    matrices = generate_matrices(len(seq_a), len(seq_b))
    filled_matrix = fill_matrix(score_mat = matrices[0], dir_mat = matrices[1], seq_a = seq_a, seq_b = seq_b, match = match, mismatch = mismatch, gap = gap)
    alns, score = find_optimal_alignment(score_mat = filled_matrix[0], dir_mat = filled_matrix[1], seq_a = seq_a, seq_b = seq_b)

    return alns, score

def generate_matrices(row: int, col: int) -> tuple[list[list[int]], list[list[list[None]]]]:
    """
    creates two (row + 1) X (col + 1) matrix with each cell as a list.
    extra row and column added to fill in initial values and corresponding direction.

    Args:
        row: length of seq_a
        col: length of seq_b

    Returns:
        generated score matrix, direction matrix
    """
    score_mat = [[ 0 for j in range(col + 1)] for i in range(row + 1)]
    dir_mat = [[ [] for j in range(col + 1)] for i in range(row + 1)]

    return score_mat, dir_mat

def fill_matrix(score_mat: list[list[int]], dir_mat:list[list[list[int]]], seq_a: str, seq_b: str, match: int, mismatch: int, gap: int) -> tuple[list[list[int]], list[list[list[int]]]]:
    """
    fill the matrix with scores based on the sequence alignment

    Args:
        score_mat: 2D initialized matrix for storing scores
        dir_mat: 2D initialized matrix for storing directions
        seq_a: first sequence
        seq_b: second sequence
        match: match score
        mismatch: mismatch penalty
        gap: gap penalty

    Returns:
        filled score matrix, direction matrix
    """
    #diagonal = 0 , right = 1, down = 2

    #fill the first row
    for col in range(1, len(seq_b) + 1):
        score_mat[0][col] = score_mat[0][col - 1] + gap
        dir_mat[0][col] = [1]

    #fill the first column
    for row in range(1, len(seq_a) + 1):
        score_mat[row][0] = score_mat[row - 1][0] + gap
        dir_mat[row][0] = [2]

    dir_mat[0][0] = [-1]

    #fill the rest of the score matrix with maximum of diagonal, right or top cell
    #fill the rest of direction matrix with the path traversed
    for row in range(1, len(seq_a) + 1):
        seq_a_base = seq_a[row - 1]
        for col in range(1, len(seq_b) + 1):
            seq_b_base = seq_b[col - 1]
            score_list = [
                score_mat[row - 1][col -1] + match if seq_a_base == seq_b_base else score_mat[row - 1][col -1] + mismatch, 
                score_mat[row][col - 1] + gap, 
                score_mat[row - 1][col] + gap
            ]
            score_mat[row][col] = max(score_list)
            dir_mat[row][col] = [i for i, j in enumerate(score_list) if j == score_mat[row][col]]

    #uncomment for printing score and direction matrices                                   
    # for row in dir_mat:
    #     print(row)    
    # for row in score_mat:
    #     print(row)
    return score_mat,dir_mat

def traverse(cell:list[int], path:int, row:int, col:int, aln_a:str, aln_b:str, dir_mat:list[list[list[int]]], seq_a:str, seq_b:str) -> tuple[str,str]:
    """
    get one optimal alignment for the sequences by backtracking the direction matrix using recursion

    Args:
        cell: current cell in direction matrix
        path: current index of cell
        row: curent row in direction matrix
        col: curent column in direction matrix
        aln_a: aligned sequence seq_a
        aln_b: aligned sequence seq_b
        dir_mat: direction matrix
        seq_a: input sequence seq_a
        seq_b: input sequence seq_b

    Returns:
        one optimal alignment
    """
    #base condition: check if cell [0,0] is reached
    if row == 0 and col == 0:
        return aln_a[::-1], aln_b[::-1]

    #recursively traverse the direction matrix to backtrack the optimal alignment
    #diagonal = 0 , right = 1, down = 2
    while path >= 0:
        #if diagonal, add base from both sequences to aln_a_temp and aln_b_temp
        if cell[path] == 0:
            return traverse(dir_mat[row-1][col-1], path, row-1, col-1, aln_a + seq_a[row-1], aln_b + seq_b[col-1], dir_mat, seq_a, seq_b)
        #if right, add base from seq_b to aln_b_temp and gap on aln_a_temp
        elif cell[path] == 1:
            return traverse(dir_mat[row][col-1], path, row, col-1, aln_a + "-", aln_b + seq_b[col-1], dir_mat, seq_a, seq_b)
        #if down, add base from seq_a to aln_a_temp and gap on aln_b_temp
        elif cell[path] == 2:
            return traverse(dir_mat[row-1][col], path, row-1, col, aln_a + seq_a[row-1], aln_b + "-", dir_mat, seq_a, seq_b)
    path -= 1


def traverse_all_paths(cell:list[int], path:int, row:int, col:int, aln_a:str, aln_b:str, dir_mat:list[list[list[int]]], seq_a:str, seq_b:str) -> tuple[str, str]:
    """
    get all optimal alignments for the sequences by backtracking the direction matrix using recursion

    Args:
        cell: current cell in direction matrix
        path: current index of cell
        row: curent row in direction matrix
        col: curent column in direction matrix
        aln_a: aligned sequence seq_a
        aln_b: aligned sequence seq_b
        dir_mat: direction matrix
        seq_a: input sequence seq_a
        seq_b: input sequence seq_b

    Returns:
        all optimal alignments
    """
    #diagonal = 0 , right = 1, down = 2
    #base condition: check if cell [0,0] is reached
    if row == 0 and col == 0:
        #print(f"{aln_a}\n{aln_b}")
        yield (aln_a[::-1], aln_b[::-1])

    #if diagonal, add base from both sequences to aln_a_temp and aln_b_temp
    if cell[path] == 0:
        new_row = row - 1
        new_col = col - 1
        cell = dir_mat[new_row][new_col]
        path = len(cell)
        aln_a_temp = aln_a + seq_a[new_row]
        aln_b_temp = aln_b + seq_b[new_col]
        #print(f"{aln_a_temp}\n{aln_b_temp}")
        while path > 0:
            yield from traverse_all_paths(cell, path - 1, new_row, new_col, aln_a_temp, aln_b_temp, dir_mat, seq_a, seq_b)
            path -= 1

    #if right, add base from seq_b to aln_b_temp and gap on aln_a_temp
    elif cell[path] == 1:
        new_col = col - 1
        cell = dir_mat[row][new_col]
        path = len(cell)
        aln_a_temp = aln_a + "-"
        aln_b_temp = aln_b + seq_b[new_col]
        #print(f"{aln_a_temp}\n{aln_b_temp}")
        while path > 0:
            yield from traverse_all_paths(cell, path - 1, row, new_col, aln_a_temp, aln_b_temp, dir_mat, seq_a, seq_b)
            path -= 1
    
    #if down, add base from seq_a to aln_a_temp and gap on aln_b_temp
    elif cell[path] == 2:
        new_row = row - 1
        cell = dir_mat[new_row][col]
        path = len(cell)
        aln_a_temp = aln_a + seq_a[new_row]
        aln_b_temp = aln_b + "-"
        #print(f"{aln_a_temp}\n{aln_b_temp}")
        while path > 0:
            yield from traverse_all_paths(cell, path - 1, new_row, col, aln_a_temp, aln_b_temp, dir_mat, seq_a, seq_b)
            path -= 1


def find_optimal_alignment(score_mat: list[list[int]], dir_mat: list[list[list[int]]], seq_a: str, seq_b: str) -> tuple[tuple[str, str], int]:
    """
    find optimal alignment for the sequences using backtracking the score matrix

    Args:
        score_mat: score matrix
        dir_mat: direction matrix
        seq_a: input sequence seq_a
        seq_b: input sequence seq_b

    Returns:
        tuple[tuple[str, str], int]: _description_
    """
    
    aln_a = ""
    aln_b = ""
    row = len(seq_a)
    col = len(seq_b)
    cell = dir_mat[row][col]
    path = len(cell) - 1

    score = score_mat[row][col]

    #to get one optimal alignment using recursion
    # a,b = traverse(cell, path, row, col, aln_a, aln_b, dir_mat, seq_a, seq_b)
    # return ((a,b), score)
    
    #to get all optimal alignments using recursion
    all_alns = [i for i in traverse_all_paths(cell, path, row, col, aln_a, aln_b, dir_mat, seq_a, seq_b)]
    
    #uncomment for printing returned alignments
    # for ele in all_alns:
    #     print(ele)

    return all_alns, score


if __name__ == "__main__":
    # alns, score = needleman_wunsch("ATAAT", "ATAATT", 1, -1, -1)
    alns, score = needleman_wunsch("TAGTCAT", "TATCAAT", 1, -1, -1)
    # alns, score = needleman_wunsch("CTTCTCGTCGGTCTCGTGGTTCGGGAAC", "CTTTCATCCACTTCGTTGCCCGGGAAC", 1, -1, -1)
    print(f"Score: {score}")
    for a in alns:
        print("\n".join(a))
        print()
    
