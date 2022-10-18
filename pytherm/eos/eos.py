from abc import ABC, abstractmethod


class EOS(ABC):
    @abstractmethod
    def get_system(self):
        pass

    @abstractmethod
    def get_f(self):
        pass

    @abstractmethod
    def get_p(self):
        pass
