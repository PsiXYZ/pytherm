from abc import abstractmethod, ABC


class activity_model(ABC):
    @abstractmethod
    def get_y(self):
        pass

    @abstractmethod
    def get_ge_nrt(self):
        pass
