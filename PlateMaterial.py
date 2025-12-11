
class Plate:
    """
    Класс пластины, характеризует её длиной - lenght, шириной - wide, толщиной - hight
    Материал - массив: [модуль юнга E, плотность материала - rho]
    """

    def __init__(self, lenght, wide, hight, material):
        # hight - толшина стержня, мала по сравнению с длинной для точности результата
        # rho - плотность материала
        self.nu = material[1] * hight * wide
        self.I = hight ** 3 * wide / 12
        self.E = material[0]
        self.lenght = lenght

    def Get_l(self):
        return self.lenght

    def Get_I(self):
        return self.I

    def Get_E(self):
        return self.E

    def Get_nu(self):
        return self.nu


class Material:
    """
    Содержит константы для материалов
    [E, rho] E - модуль юнга, rho - плотность
    """

    def __init__(self, E, rho):
        self.mat = [E, rho]

    def Get_Material(self):
        return self.mat

    @classmethod
    def Steel(cls):  # Сталь
        # 0 - Модуль юнга E в паскалях, второе - rho: плотность в кг / м^3
        return [200 * 10 ** 9, 8000]

    @classmethod
    def Titan(cls):  # Титан
        # 0 - Модуль юнга E в паскалях, второе - rho: плотность в кг / м^3
        return [120 * 10 ** 9, 7870]

    @classmethod
    def Aluminum(cls):  # Алюминий
        # 0 - Модуль юнга E в паскалях, второе - rho: плотность в кг / м^3
        return [70 * 10 ** 9, 2700]

    @classmethod
    def Copper(cls):  # Медь
        # 0 - Модуль юнга E в паскалях, второе - rho: плотность в кг / м^3
        return [125 * 10 ** 9, 8920]

    @classmethod
    def Brass(cls):  # латунь
        # 0 - Модуль юнга E в паскалях, второе - rho: плотность в кг / м^3
        return [95 * 10 ** 9, 8500]
