
import pandas as pd
import numpy as np


metadata_path_lotus = '/Users/pma/210523_frozen_metadata.csv'

metadata_path_dnp = '/Users/pma/210523_dnp_metadata.csv'


dt_metadata_lotus = pd.read_csv(metadata_path_lotus,
                                sep=',', error_bad_lines=False, low_memory=True)


dt_metadata_dnp = pd.read_csv(metadata_path_dnp,
                                sep=',', error_bad_lines=False, low_memory=True)



dt_metadata = pd.concat([dt_metadata_lotus, dt_metadata_dnp], ignore_index=True, sort=False)

# Since we have structure-organisms pairs duplicated across references in the lotus set we drop them
dt_metadata_deduped = dt_metadata.drop_duplicates(['structure_inchikey', 'organism_name'])

dt_metadata_lotus.info()
dt_metadata_dnp.info()
dt_metadata.info()
dt_metadata_deduped.info()

dt_metadata_deduped.to_csv('/Users/pma/210523_lotus_dnp_metadata.csv', index = None)