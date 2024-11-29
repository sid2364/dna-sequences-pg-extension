DROP TABLE IF EXISTS kmer_data_t;
DO $$
BEGIN
    IF NOT EXISTS (SELECT 1 FROM information_schema.tables WHERE table_name = 'kmer_data_t') THEN
        CREATE TABLE kmer_data_t (
            id SERIAL PRIMARY KEY,
            kmer_a kmer
        );
    END IF;
END $$;


DO $$
BEGIN
    FOR i IN 1..100 LOOP
        INSERT INTO kmer_data_t (kmer_a) VALUES ('ATC');
    END LOOP;
END $$;

EXPLAIN ANALYZE
SELECT * FROM kmer_data_t WHERE kmer_a = 'ATC';


DROP TABLE IF EXISTS kmer_data_t;
DO $$
BEGIN
    IF NOT EXISTS (SELECT 1 FROM information_schema.tables WHERE table_name = 'kmer_data_t') THEN
        CREATE TABLE kmer_data_t (
            id SERIAL PRIMARY KEY,
            kmer_a kmer
        );
    END IF;
END $$;

CREATE INDEX spgist_kmer_idx
ON kmer_data_t USING spgist (kmer_a spgist_kmer_ops);

SET enable_seqscan = OFF;
EXPLAIN ANALYZE
SELECT * FROM kmer_data_t WHERE kmer_a = 'ATC';