import pandas as pd
import pyopenms

db = pyopenms.ModificationsDB()

mods = []

db.getAllSearchModifications(mods)

data = []

for mod in mods:
    try:
        row = {'modification': mod}

        idx = db.findModificationIndex(mod)

        res_mod = db.getModification(idx)

        row['weight'] = res_mod.getDiffFormula().getMonoWeight()

        row['terminus_specificity'] = res_mod.getTermSpecificityName(res_mod.getTermSpecificity())

        data.append(row)

    except RuntimeError:
        pass

data = pd.DataFrame(data, columns=['modification', 'weight', 'terminus_specificity'])

data.to_csv('../soil/data/protein_mods.tsv', index=False, sep='\t')
