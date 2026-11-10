"""Read penultimate LAST identity; independently validate final alignment parameters."""
import re
from pathlib import Path

PREFIX = '# substitution percent identity:'
IDENTITY_SOURCE = 'last-train-penultimate'


def parse_identity_values(path):
    lines = Path(path).read_text().splitlines()
    indices = [i for i, line in enumerate(lines) if line.startswith(PREFIX)]
    if len(indices) < 2:
        raise ValueError(f'At least two LAST identity markers required (penultimate identity): {path}')
    i = indices[-1]
    values = []
    for position, index in (('penultimate', indices[-2]), ('final', i)):
        text = lines[index][len(PREFIX):].strip()
        if not re.fullmatch(r'[0-9]+(?:[.][0-9]+)?', text) or not 0 <= float(text) <= 100:
            raise ValueError(f'Invalid {position} LAST identity: {path}')
        values.append(float(text))
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
    return tuple(values)


def parse_identity(path):
    return parse_identity_values(path)[0]


if __name__ == '__main__':
    import sys
    print(IDENTITY_SOURCE if sys.argv[1:] == ['--source'] else parse_identity(sys.argv[1]))
