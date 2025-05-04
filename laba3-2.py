import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm

sns.set_theme(style="ticks")


class TimeSeriesPlotter:
    def __init__(self, dataset_name: str):  # Исправлено на __init__
        self.dataset_name = dataset_name
        self.data = self.load_data()

    def load_data(self):
        """Загружает и возвращает DataFrame из statsmodels"""
        dataset = getattr(sm.datasets, self.dataset_name).load_pandas()
        df = dataset.data.copy()
        df.index.name = "date"  # Даём имя индексу для удобства
        df = df.reset_index()  # Переносим дату в отдельную колонку
        return df

    def prepare_data(self, date_column: str, value_column: str):
        """Приводит данные к формату: дата + значения"""
        if date_column not in self.data.columns or value_column not in self.data.columns:
            raise ValueError(f"Столбец {date_column} или {value_column} отсутствует в данных.")
        
        self.data[date_column] = pd.to_datetime(self.data[date_column])
        df = self.data[[date_column, value_column]].dropna()
        self.data = df.rename(columns={date_column: "date", value_column: "value"})

    def filter_date_range(self, start: str, end: str):
        """Фильтрация данных по временным рамкам"""
        self.data = self.data[(self.data["date"] >= start) & (self.data["date"] <= end)]

    def plot(self, title: str = "Временной ряд", ylabel: str = "Значение"):
        """Строит временной график с помощью Seaborn"""
        plt.figure(figsize=(12, 6))
        ax = sns.lineplot(
            data=self.data,
            x="date",
            y="value",
            label=self.dataset_name,
            color="darkblue",
        )
        ax.set_title(title, fontsize=14, fontweight="bold")
        ax.set_xlabel("Дата", fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.grid(True)
        plt.xticks(rotation=30)  # наклонные даты
        plt.legend()
        plt.tight_layout()  # гармоничное расположение
        plt.show()
# Создаём экземпляр класса с корректным названием датасета
plotter = TimeSeriesPlotter("co2")

# Проверяем загруженные данные
print(plotter.data.head())

# Преобразуем данные для построения графика
plotter.prepare_data(date_column="date", value_column="co2")

# Фильтруем временной диапазон
plotter.filter_date_range("1958-01-01", "1980-12-31")

# Строим график
plotter.plot(title="Динамика уровня CO2 (1958–1980)", ylabel="CO2 (ppm)")
