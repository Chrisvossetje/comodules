import itertools

def compose(p1, p2):
    """Return the composition of two permutations: p1 after p2"""
    return tuple(p1[i] for i in p2)

def permutation_multiplication_table(n):
    perms = list(itertools.permutations(range(n)))
    index_map = {p: i for i, p in enumerate(perms)}
    
    size = len(perms)
    table = [[0]*size for _ in range(size)]
    
    for i, p1 in enumerate(perms):
        for j, p2 in enumerate(perms):
            composed = compose(p1, p2)
            table[i][j] = index_map[composed]
    
    return perms, table

def print_table(perms, table):
    print("Permutation Multiplication Table (entries are indices of the resulting permutations):\n")
    print("       ", end="")
    for i in range(len(perms) - 2):
        print(f"   {chr(i + 97 + 2)}", end="")
    print()
    for i, row1 in enumerate(table):
        for j, row2 in enumerate(table):
            if (j <= 1):
                continue
            print(f"{(chr(i+97)):2} * {chr(j+97):2}", end="")
            for k, _ in enumerate(table):
                if k <= 1:
                    continue
                if k == table[i][j]:
                    print(f"   1", end="")
                else:
                    print(f"   0", end="")

            print()

def perm_to_string(perm):
    out = ""
    for p in perm:
        out += str(p)
    return out

if __name__ == "__main__":
    n = 3  # Change this for other sizes
    perms, table = permutation_multiplication_table(n)
    
    print("Permutations:")
    for i, p in enumerate(perms):
        print(f"{i}: {perm_to_string(p)}")
    
    print()
    print_table(perms, table)