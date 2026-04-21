import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# === 1. Считывание файла ===
filename = r"C:\Users\6anna\CLionProjects\Sem6Laba3_FDTD\cmake-build-debug\field_Gauss_soft_PML.csv"
df = pd.read_csv(filename)

print(df.head())
print(df.columns)

# Ожидаются колонки:
# time_over_fL, x_tilde, Ez

# === 2. Преобразование в 2D-матрицу ===
# строки -> время, столбцы -> координата x
pivot = df.pivot(index="time_over_fL", columns="x_tilde", values="Ez")

times = pivot.index.to_numpy()
x = pivot.columns.to_numpy()
Ez = pivot.to_numpy()

# === 3. Heatmap поля Ez(x, t) ===
plt.figure(figsize=(10, 6))
plt.imshow(
    Ez,
    aspect='auto',
    origin='lower',
    extent=[x.min(), x.max(), times.min(), times.max()],
    cmap='RdBu_r'
)
plt.colorbar(label='Ez')
plt.xlabel('x_tilde')
plt.ylabel('time_over_fL')
plt.title('Поле Ez(x, t)')
plt.tight_layout()
plt.show()

# === 4. Несколько пространственных срезов Ez(x) ===
plt.figure(figsize=(10, 6))

n_slices = 6
indices = np.linspace(0, len(times) - 1, n_slices, dtype=int)

for idx in indices:
    plt.plot(x, Ez[idx], label=f"t = {times[idx]:.3f}")

plt.xlabel('x_tilde')
plt.ylabel('Ez')
plt.title('Пространственные срезы поля Ez(x)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# === 5. Временной срез в одной точке x (например, около центра) ===
x_index = len(x) // 2

plt.figure(figsize=(10, 5))
plt.plot(times, Ez[:, x_index])
plt.xlabel('time_over_fL')
plt.ylabel('Ez')
plt.title(f'Временная зависимость Ez(t) в точке x = {x[x_index]:.3f}')
plt.grid(True)
plt.tight_layout()
plt.show()