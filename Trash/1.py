from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow

import sys


class Window(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("sdasdasd")
        self.setGeometry(300, 350, 350, 200)

        self.new_text = QtWidgets.QLabel(self)



        self.main_text = QtWidgets.QLabel(self)
        self.main_text.setText('Это ндапись')
        self.main_text.move(100, 100)
        self.main_text.adjustSize()

        self.btn = QtWidgets.QPushButton(self)
        self.btn.move(70, 150)
        self.btn.setText('Нажми')
        self.btn.setFixedWidth(200)
        self.btn.clicked.connect(self.add_label)


    def add_label(self):
        self.new_text.setText('Вторая надпись')
        self.new_text.move(50,50)
        self.new_text.adjustSize()

def aplication():
    app = QApplication(sys.argv)
    window = Window()

    window.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    aplication()