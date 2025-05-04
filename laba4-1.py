from Bio import SeqIO


# Класс для обработки GenBank-файлов
class GenBankProcessor:
    def init(self, input_path, output_path):
        self.input_path = input_path  # Путь к исходному GenBank-файлу (sequence.gb)
        self.output_path = output_path  # Путь для сохранения объединённого файла
        self.records = []  # Список для хранения всех считанных записей

    def load_records(self):
        try:
            # Считываем все записи и преобразуем их в список
            self.records = list(SeqIO.parse(self.input_path, "genbank"))
            print(f"Загружено {len(self.records)} записей.")
        except FileNotFoundError:
            # Обрабатываем ошибку, если файл не найден
            print(f"Файл {self.input_path} не найден.")
            raise  # Повторно выбрасываем исключение, чтобы программа завершилась

    def count_cds(self):
        cds_count = 0
        for record in self.records:
            # Проходим по всем признакам (features) записи и считаем те, что типа 'CDS'
            cds_count += sum(1 for feature in record.features if feature.type == "CDS")
        return cds_count  # Возвращаем общее количество CDS

    def save_records(self):
        # Сохраняем загруженные записи в новый GenBank-файл
        with open(self.output_path, "w") as out_handle:

            SeqIO.write(self.records, out_handle, "genbank")
        print(f"Сохранено в файл: {self.output_path}")


# Создаём объект GenBankProcessor, указывая входной и выходной файл
processor = GenBankProcessor("sequence.gb", "combined_file.gbk")

# Загружаем записи из файла
processor.load_records()

# Считаем количество CDS и выводим результат
cds_total = processor.count_cds()
print(f"Общее количество CDS: {cds_total}")

# Сохраняем все записи в новый файл
processor.save_records()