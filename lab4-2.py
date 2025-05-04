from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.Seq import UndefinedSequenceError


def calculate_gc_content(record):
    try:
        if not record.seq:  
            return 0.0
        return gc_fraction(record.seq)
    except UndefinedSequenceError:
        return 0.0 


def analyze_gc_content(input_file):
    try:
        records = list(SeqIO.parse(input_file, "genbank"))
        if not records:
            print("Файл не содержит записей GenBank!")
            return

        gc_data = []
        for record in records:
            try:
                gc_content = calculate_gc_content(record)
                description = record.description
                record_id = record.id
                gc_data.append((gc_content, description, record_id))
            except Exception as e:
                print(f"Ошибка при обработке записи {record.id}: {str(e)}")
                continue

        if not gc_data:
            print("Не удалось обработать ни одной записи!")
            return

        gc_data.sort()  

        for gc, desc, rec_id in gc_data:
            print(f"{rec_id}: {desc}, GC = {gc:.4f}")  

    except FileNotFoundError:
        print(f"Файл {input_file} не найден!")
    except Exception as e:
        print(f"Критическая ошибка: {str(e)}")


if __name__ == "__main__":
    input_file = "merged_sequences.gb"  
    analyze_gc_content(input_file)