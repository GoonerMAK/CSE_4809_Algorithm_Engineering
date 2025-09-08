# Correct 2D Rabin-Karp implementation and test with the user's example
def rabin_karp_2d(text, pattern, base_col=31, base_row=37, mod=10**9+7):
    """
    2D Rabin-Karp search:
    - text: list of strings (R rows, each length C)
    - pattern: list of strings (r rows, each length c)
    Returns list of top-left matches (i, j) using 0-based indexing.
    """
    R = len(text)
    C = len(text[0]) if R > 0 else 0
    r = len(pattern)
    c = len(pattern[0]) if r > 0 else 0

    if r == 0 or c == 0 or R < r or C < c:
        return []

    # helper to convert char to integer value (consistent across text & pattern)
    def val(ch):
        return ord(ch)  # using raw ord ensures we handle arbitrary characters

    # Precompute powers
    pow_col = pow(base_col, c - 1, mod)
    pow_row = pow(base_row, r - 1, mod)

    # 1) Compute pattern hash
    pattern_row_hashes = []
    for row in pattern:
        h = 0
        for ch in row:
            h = (h * base_col + val(ch)) % mod
        pattern_row_hashes.append(h)

    pattern_hash = 0
    for rh in pattern_row_hashes:
        pattern_hash = (pattern_hash * base_row + rh) % mod

    # 2) Compute rolling row hashes for text (for every row, for every column window)
    row_hashes = [[0] * (C - c + 1) for _ in range(R)]
    for i in range(R):
        # initial hash for row i
        h = 0
        for j in range(c):
            h = (h * base_col + val(text[i][j])) % mod
        row_hashes[i][0] = h
        # rolling across columns
        for j in range(1, C - c + 1):
            left_val = val(text[i][j - 1])
            right_val = val(text[i][j + c - 1])
            h = (h - left_val * pow_col) % mod
            h = (h * base_col + right_val) % mod
            row_hashes[i][j] = h

    matches = []

    # 3) For each column window j, compute vertical rolling hashes using row_hashes[*][j]
    for j in range(C - c + 1):
        # initial vertical hash for first r rows at column j
        h = 0
        for i in range(r):
            h = (h * base_row + row_hashes[i][j]) % mod

        # check and slide down
        if h == pattern_hash:
            # verify to avoid false positives from hash collisions
            match = True
            for di in range(r):
                if text[di][j:j+c] != pattern[di]:
                    match = False
                    break
            if match:
                matches.append((0, j))

        for i in range(1, R - r + 1):
            # remove top row contribution and add new bottom row
            top_hash = row_hashes[i - 1][j]
            bottom_hash = row_hashes[i + r - 1][j]
            h = (h - top_hash * pow_row) % mod
            h = (h * base_row + bottom_hash) % mod

            if h == pattern_hash:
                # verify to avoid false positives from hash collisions
                match = True
                for di in range(r):
                    if text[i + di][j:j + c] != pattern[di]:
                        match = False
                        break
                if match:
                    matches.append((i, j))

    return matches

# User's provided example
text = [
    "abcdabc",
    "bcpikbc",
    "cduuucd",
    "daapika",
    "aabuuub"
]

pattern = [
    "pik",
    "uuu"
]

matches = rabin_karp_2d(text, pattern)
print("Matches found at:", matches)

# For demonstration, print the matched region if any
if not matches:
    print("No match found")
else:
    for (i, j) in matches:
        print(f"Match at top-left = ({i}, {j}):")
        for di in range(len(pattern)):
            print(text[i + di][j:j + len(pattern[0])])
