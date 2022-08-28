
import seaborn as sns

from collections import defaultdict

def get_nth(palette, n):
    return sns.color_palette(palette)[n]

BATCH_PALETTE = defaultdict(lambda : 'lightgrey', {
    's4d1' : get_nth('light:blue', 1),
    's4d9' : get_nth('light:blue', 2),
    's4d8' : get_nth('light:blue', 3),
    's2d1' : get_nth('light:green', 1),
    's2d4' : get_nth('light:green', 2),
    's2d5' : get_nth('light:green', 3),
    's3d3' : get_nth('light:red', 4),
    's3d6' : get_nth('light:red', 3),
    's3d7' : get_nth('light:red', 1),
    's3d10' : get_nth('light:red', 2),
    's1d1' : get_nth('light:orange', 1),
    's1d2' : get_nth('light:orange', 2),
    's1d3' : get_nth('light:orange', 3),
})

CELL_PALETTE = {'B1 B': '#023fa5',
 'CD4+ T activated': '#7d87b9',
 'CD4+ T naive': '#d6bcc0',
 'CD8+ T': '#bec1d4',
 'CD8+ T naive': '#bb7784',
 'CD14+ Mono': '#8e063b',
 'CD16+ Mono': '#4a6fe3',
 'Erythroblast': '#8595e1',
 'G/M prog': '#b5bbe3',
 'HSC': '#e6afb9',
 'ID2-hi myeloid prog': '#e07b91',
 'ILC': '#d33f6a',
 'Lymph prog': '#11c638',
 'MK/E prog': '#8dd593',
 'NK': '#c6dec7',
 'Naive CD20+ B': '#ead3c6',
 'Normoblast': '#f0b98d',
 'Plasma cell': '#ef9708',
 'Proerythroblast': '#0fcfc0',
 'Transitional B': '#9cded6',
 'cDC2': '#d5eae7',
 'pDC': '#f3e1eb'}