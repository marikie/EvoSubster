"""Parse only LAST's complete final training block."""
import re
from pathlib import Path

PREFIX = '# substitution percent identity:'


def parse_identity(path):
    lines = Path(path).read_text().splitlines()
    indices = [i for i, line in enumerate(lines) if line.startswith(PREFIX)]
    if not indices:
        raise ValueError(f'Missing LAST identity: {path}')
    i = indices[-1]
    text = lines[i][len(PREFIX):].strip()
    if not re.fullmatch(r'[0-9]+(?:[.][0-9]+)?', text) or not 0 <= float(text) <= 100:
        raise ValueError(f'Invalid LAST identity: {path}')
    block = lines[i+1:]
    if not all(any(line.startswith('#last -' + option) for line in block)
               for option in ('t', 'a', 'A', 'b', 'B', 'S')):
        raise ValueError(f'Incomplete LAST parameters: {path}')
    header = '# score matrix (query letters = columns, reference letters = rows):'
    if header not in block:
        raise ValueError(f'Missing LAST score matrix: {path}')
    rows = block[len(block) - 1 - block[::-1].index(header) + 1:]
    number = r'-?[0-9]+(?:[.][0-9]+)?'
    pattern = r'^([ACGT])\s+' + r'\s+'.join([number] * 4) + r'\s*$'
    letters = [match.group(1) for line in rows if (match := re.match(pattern, line))]
    if set(letters) != set('ACGT') or len(letters) != 4:
        raise ValueError(f'Incomplete LAST score matrix: {path}')
    return float(text)


if __name__ == '__main__':
    import sys
    print(parse_identity(sys.argv[1]))
