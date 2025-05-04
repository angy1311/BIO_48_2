from Bio import SeqIO


def translate_cds_features(genbank_file):
    for record in SeqIO.parse(genbank_file, "genbank"):
        print(f"\n{record.id}: {record.description}")

        cds_features = [feature for feature in record.features if feature.type == "CDS"]

        if not cds_features:
            print("Нет кодирующих областей (CDS) в этой записи")
            continue

        for feature in cds_features:
            location = feature.location
            strand = "+" if location.strand > 0 else "-"

            print(f"Coding sequence location = [{location.start}:{location.end}]({strand})")

            try:
                cds_seq = location.extract(record.seq)
                protein_seq = cds_seq.translate(to_stop=True)

                print("Translation =")
                for i in range(0, len(protein_seq), 60):
                    print("  ", str(protein_seq[i:i + 60]))

            except Exception as e:
                print(f"Ошибка при трансляции: {str(e)}")
                continue


if __name__ == "__main__":
    input_file = "merged_sequences.gb"  
    translate_cds_features(input_file)