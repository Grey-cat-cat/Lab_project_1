import matplotlib.animation
import matplotlib.pyplot
import numpy as np
import Equations_finite


class Gif_Animation:
    """
    equation содержит данные о колебании пластины во времени,
    они используются для создания анимации
    Класс создаёт анимацию и сохраняет её
    Для 
    """

    def __init__(self, equation: Equations_finite.Euler_Bernoulli_Dinamic_Finite, nubber_of_frames, gif_name='plate_animation', fps=10):
        """
        equation - объект класса Euler_Bernoulli_Dinamic, nubber_of_frames - количество кадров в итоговой анимации
        gif_name - название Gif файла анимации, fps - количество кадров в секунду
        """
        matplotlib.pyplot.close('all')
        matplotlib.pyplot.pause(0.1)

        self.gif_name = gif_name
        self.eq = equation
        self.frame_n = nubber_of_frames
        self.delay_milsec = 1000 / fps
        self.fps = fps

        # Создание объекта, в котором создаётся график, и осей
        self.figure, self.axis = matplotlib.pyplot.subplots(figsize=(12, 12))

        # Максимальное значение в графиках
        all_data = self.eq.result[:self.frame_n, :self.eq.N_x]
        self.y_max = np.max(all_data) * 1.1
        self.y_min = -self.y_max
        self.line = None

    def Init(self):
        """Начало анимации, вызывается 1 раз"""
        self.axis.clear()

        self.axis.set_xlim(0, self.eq.L)
        self.axis.set_ylim(self.y_min, self.y_max)

        self.axis.set_title('Колебания')
        self.axis.set_xlabel('Координата точки, м')
        self.axis.set_ylabel('Прогиб, мм')

        self.line, = self.axis.plot(
            [], [], 'b-', linewidth=2, label='m - x; mm - y')

        return (self.line,)

    def Update(self, frame):
        """
        Покадрово обновляет график, создавая анимацию.
        frame - номер времени в симуляции у Euler_Bernoulli_Dinamic
        """
        if frame >= self.frame_n:
            return (self.line,)

        array_x = np.zeros(self.eq.N_x)
        array_y = np.zeros(self.eq.N_x)

        array_y = self.eq.result[frame][:self.eq.N_x]

        for i in range(self.eq.N_x):  # Координата x
            array_x[i] = self.eq.dx * i

        # Время
        current_time = frame * self.eq.T / (self.eq.N_t)
        self.axis.set_title(f'Время: {current_time:.3f} с')

        self.line.set_data(array_x, array_y)

        return (self.line,)

    def save_gif(self):
        """Создаёт и сохраняет анимацию в виде Gif файла"""
        print("started creating Gif")

        animation = matplotlib.animation.FuncAnimation(
            fig=self.figure, func=self.Update, frames=self.frame_n, init_func=self.Init, interval=self.delay_milsec, repeat=True)

        animation.save(self.gif_name, writer='pillow', fps=self.fps, dpi=100)
        print(f"Gif сохранён под именем: {self.gif_name}")
        matplotlib.pyplot.show()
        matplotlib.pyplot.close(self.figure)
