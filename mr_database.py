import sqlite3

def store_instruments_to_db(instruments):
    conn = sqlite3.connect('instruments.db')
    cursor = conn.cursor()
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS instruments (
            chr_exposure TEXT,
            se_exposure REAL,
            pval_exposure REAL,
            beta_exposure REAL,
            pos_exposure INTEGER,
            samplesize_exposure INTEGER,
            id_exposure TEXT,
            SNP TEXT,
            effect_allele_exposure TEXT,
            other_allele_exposure TEXT,
            eaf_exposure REAL,
            exposure TEXT,
            mr_keep_exposure INTEGER,
            pval_origin_exposure REAL,
            data_source_exposure TEXT
        )
    ''')
    for instrument in instruments:
        cursor.execute('''
            INSERT INTO instruments (
                chr_exposure, se_exposure, pval_exposure, beta_exposure, pos_exposure, samplesize_exposure,
                id_exposure, SNP, effect_allele_exposure, other_allele_exposure, eaf_exposure, exposure,
                mr_keep_exposure, pval_origin_exposure, data_source_exposure
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            instrument['chr.exposure'], instrument['se.exposure'], instrument['pval.exposure'], instrument['beta.exposure'],
            instrument['pos.exposure'], instrument['samplesize.exposure'], instrument['id.exposure'], instrument['SNP'],
            instrument['effect_allele.exposure'], instrument['other_allele.exposure'], instrument['eaf.exposure'],
            instrument['exposure'], instrument['mr_keep.exposure'], instrument['pval_origin.exposure'], instrument['data_source.exposure']
        ))
    conn.commit()
    conn.close()

def get_all_instruments():
    conn = sqlite3.connect('instruments.db')
    cursor = conn.cursor()
    cursor.execute('SELECT * FROM instruments')
    columns = [column[0] for column in cursor.description]
    rows = cursor.fetchall()
    conn.close()
    instruments = [dict(zip(columns, row)) for row in rows]
    return instruments

def get_snps_from_db():
    conn = sqlite3.connect('instruments.db')
    cursor = conn.cursor()
    cursor.execute('SELECT SNP FROM instruments')
    snps = [row[0] for row in cursor.fetchall()]
    conn.close()
    return snps