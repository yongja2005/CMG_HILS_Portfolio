import sys
import time

from PyQt5.QtWidgets import *
from PyQt5 import uic
from PyQt5.QtCore import *
from PyQt5.QtGui import QPainter
from PyQt5.QtChart import QLineSeries, QChart, QValueAxis, QDateTimeAxis
from PyQt5.QtCore import Qt, QDateTime
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import QFile
import serial
import numpy as np
import pandas as pd
import math
import Utills

port_number = 'COM6'
baud_rate = '115200'
ser = serial.Serial( port_number , baud_rate)
form_class = uic.loadUiType("cmg_hils.ui")[0]

class CMG(QThread):
  # pyqtSignal(type)에서 type을 @pyqtSlot(type)과 꼭 일치 시켜줘야함
  finished = pyqtSignal(np.ndarray)
  finished2 = pyqtSignal(np.ndarray)
  finished3 = pyqtSignal(np.ndarray)

  def run(self):
    # stepping motor conf.
    step_angle = Utills.degtorad(0.05625)  # min. step per 0.1sec

    # Initial condition
    dt = 0.1  # sec(delta t)
    time = np.arange(0, 50 + dt, dt)
    I_sc = np.array([
      [0.028, 0, 0], [0, 0.028, 0], [0, 0, 0.034]
    ], float)  # kgm^2
    beta = Utills.degtorad(53.13)  # deg to rad
    h_a = 1.312e-6 * (100 * math.pi)  # Nmsec (100pi rad/s = 3000 RPM)
    zeta = 0.707  # damping ratio
    t_s = 15  # settling time(sec)
    w_n = 4 / zeta / t_s  # natural frequency
    k_D = 2 * zeta * w_n * I_sc  # D gain
    k_P = 2 * w_n * w_n * I_sc  # P gain
    limit_GV = math.pi  # rad/s
    e_0 = 0.01

    # create empty array
    w = np.zeros(shape=(len(time), 3))
    q_cur = np.zeros(shape=(len(time), 4))
    u_cmd = np.zeros(shape=(len(time), 3))
    gamma = np.zeros(shape=(len(time), 4))
    dotgamma = np.zeros(shape=(len(time), 4))
    dotgamma_real = np.zeros(shape=(len(time), 4))
    h_cmg = np.zeros(shape=(len(time), 3))
    u_cmg = np.zeros(shape=(len(time), 3))
    u_real = np.zeros(shape=(len(time), 3))
    RPY = np.zeros(shape=(len(time), 3))
    command_step = np.zeros(shape=(len(time), 4))
    m = np.zeros(shape=(len(time), 1))

    w[0, :] = Utills.degtorad(0, 0, 0)  # deg/s to rad/s (max = 0.1rad/s)
    euler_0 = Utills.degtorad(0, 0, 0)  # deg to rad
    euler_cmd = Utills.degtorad(-36, 20, 30)
    q_cur[0, :] = Utills.transf(euler_0[0, :])
    q_cmd = Utills.transf(euler_cmd[0, :])
    gamma[0, :] = Utills.degtorad(0, 0, 0, 0)  # deg to rad
    RPY[0, :] = Utills.transf(q_cur[0, :])

    while True:
      for k in range(len(time) - 1):
        # Command torque
        q_e = Utills.skew(q_cmd[0, :]).dot(q_cur[k, :])
        u_cmd[k, :] = np.array([
          - (k_P.dot(q_e[0:3].T) + k_D.dot(w[k, :].T))
        ], float)

        # Angular Momentum about 3-axis of cmg : h_cmg
        h_cmg[k, :] = Utills.AngMomCMG(gamma[k, :], beta, h_a)

        # Torque about 3-axis of cmg : u_cmg
        u_cmg[k, :] = np.array([
          -(u_cmd[k, :].T + np.cross(w[k, :], h_cmg[k, :]))
        ], float)

        # Gimbal velocity
        tt = dt * k
        E = Utills.dither(tt, e_0)
        dotA = Utills.dotconfMAT(gamma[k, :], beta)
        m[k,0] = np.linalg.det(dotA.dot(dotA.T))
        a = e_0 * np.exp(-10 * m[k,0])
        dotAA = dotA.T.dot(np.linalg.inv(dotA.dot(dotA.T) + a * E))
        dotgamma[k, :] = np.array([
          dotAA.dot(u_cmg[k, :].T) / h_a
        ])

        for i in range(0, 4):  # dotgamma를 step(0.1125deg)별로 쪼개기
          if dotgamma[k, i] > limit_GV:
            temp = (limit_GV * dt) // step_angle
            dotgamma_real[k, i] = temp * step_angle / dt
          elif dotgamma[k, i] < -limit_GV:
            temp = -(limit_GV * dt) // step_angle
            dotgamma_real[k, i] = temp * step_angle / dt
          elif -step_angle < dotgamma[k, i] * dt < step_angle:
            dotgamma_real[k, i] = 0
          else:
            if dotgamma[k, i] > 0:
              temp = dotgamma[k, i] * dt // step_angle
              dotgamma_real[k, i] = temp * step_angle / dt
            else:
              temp = dotgamma[k, i] * dt // step_angle
              dotgamma_real[k, i] = (temp + 1) * step_angle / dt

        gamma[k + 1, :] = np.array([gamma[k, :] + dotgamma_real[k, :] * dt], float)
        command_step[k, :] = np.round(dotgamma_real[k, :] / step_angle * dt)

        Trans = str(command_step[k, 0]) + ',' + str(command_step[k, 0]) + 'Q'
        Trans = Trans.encode('utf-8')
        ser.write(Trans)
        self.msleep(100)

        # feedback from encoder
        gamma[k+1,0] = ser.readline().decode("utf-8")
        dotgamma_real[k,:] = np.array([ ( gamma[k+1,:] - gamma[k,:] )/dt ], float)

        self.finished.emit(RPY[k, :])
        self.finished2.emit(dotgamma_real[k, :])
        self.finished3.emit(gamma[k, :])

        # Real torque of cmg : u_real
        u_real[k, :] = np.array([
          -(h_a * dotA.dot(dotgamma_real[k, :]) + np.cross(w[k, :], h_cmg[k, :]))
        ], float)

        # Dynamic Eqauation
        w[k + 1, :] = Utills.w_rk45(w[k, :], u_real[k, :], -1 + k, dt, I_sc)

        # Kinematics Equation
        q_cur[k + 1, :] = Utills.q_rk45(w[k + 0, :], q_cur[k + 0, :], -1 + k, dt)

        # Quaternion -> euler angle(deg)
        RPY[k + 1, :] = Utills.transf(q_cur[k + 1, :])

      print('finish')
      # pandas DataFrame으로 변경 후 csv확장자로 저장
      RPY_close = pd.DataFrame(RPY)
      RPY_close.to_csv('RPY_close.csv', index=False, header=None)
      w_close = pd.DataFrame(w)
      w_close.to_csv('w_close.csv', index=False, header=None)
      dotgamma_real_close = pd.DataFrame(dotgamma_real)
      dotgamma_real_close.to_csv('dotgamma_real_close.csv', index=False, header=None)
      gamma_close = pd.DataFrame(gamma)
      gamma_close.to_csv('gamma_close.csv', index=False, header=None)
      u_real_close = pd.DataFrame(u_real)
      u_real_close.to_csv('u_real_close.csv', index=False, header=None)
      m_close = pd.DataFrame(m)
      m_close.to_csv('m_close.csv', index=False, header=None)
      break

class MyWindow(QMainWindow, form_class, QWidget):
  def __init__(self):
    super().__init__()
    self.setupUi(self)
    self.setWindowTitle('ODD\'s CMG HILS')
    self.setWindowIcon(QIcon("favicon.png"))

    self.viewLimit = 500
    self.last_viewLimit = int(self.viewLimit / 10 + 10)

    self.Roll = QLineSeries()
    self.Pitch = QLineSeries()
    self.Yaw = QLineSeries()

    self.EulerChart = QChart()
    self.EulerChart.addSeries(self.Roll)
    self.EulerChart.addSeries(self.Pitch)
    self.EulerChart.addSeries(self.Yaw)

    axisX = QDateTimeAxis()
    axisX.setFormat("HH:mm:ss")
    axisX.setTickCount(7)
    dt = QDateTime.currentDateTime()
    axisX.setRange(dt, dt.addSecs(self.viewLimit))
    axisY = QValueAxis()

    # 차트 범례
    self.EulerChart.legend().markers(self.Roll)[0].setLabel("roll")
    self.EulerChart.legend().markers(self.Pitch)[0].setLabel("pitch")
    self.EulerChart.legend().markers(self.Yaw)[0].setLabel("yaw")
    # 축 레이블
    axisY.setVisible()

    self.EulerChart.addAxis(axisX, Qt.AlignBottom)
    self.EulerChart.addAxis(axisY, Qt.AlignLeft)
    self.Roll.attachAxis(axisX)
    self.Roll.attachAxis(axisY)
    self.Pitch.attachAxis(axisX)
    self.Pitch.attachAxis(axisY)
    self.Yaw.attachAxis(axisX)
    self.Yaw.attachAxis(axisY)

    self.EulerChart.layout().setContentsMargins(0, 0, 0, 0)

    self.Eulerangle.setChart(self.EulerChart)
    self.Eulerangle.setRenderHints(QPainter.Antialiasing)

    self.label_6.setText(port_number)
    self.label_8.setText(baud_rate)

    self.cmg = CMG()
    self.pushButton_2.clicked.connect(self.btn_clicked_2)

  @pyqtSlot(np.ndarray)
  def Current_RPY(self, RPY):
    try:
      self.lineEdit.setText(str(round(RPY[2], 2)))  # Roll
      self.lineEdit_2.setText(str(round(RPY[1], 2)))  # Pitch
      self.lineEdit_3.setText(str(round(RPY[0], 2)))  # Yaw

      # 현재 시간 정보를 얻어와서 시간과 현재가를 함께 저장. append 메서드는 millisecond를 입력받기 때문에 MSecsSinceEpoch() 메서드로 QDateTime 객체를 millisecond로 변환.
      dt = QDateTime.currentDateTime()
      self.Roll.append(dt.toMSecsSinceEpoch(), RPY[0])
      self.Pitch.append(dt.toMSecsSinceEpoch(), RPY[1])
      self.Yaw.append(dt.toMSecsSinceEpoch(), RPY[2])
      # 차트의 축정보를 업데이트하는 _updateAxis() 메서드를 호출. 실시간으로 추가되는 데이터의 위치를 지정.
      self.__updateAxis()
    except:
      pass

  def __updateAxis(self):
    try:
      # pointsVector 메서드를 사용해서 QLineSeries 객체에 저장된 데이터를 리스트로 얻어옴. pvs에 저장된 리스트 안에는 QPointF 객체로 위치 정보가 저장되어 있음
      pvs = self.Roll.pointsVector()
      # 꺼내온 x좌표 데이터가 ms라서 fromMSecsSinceEpoch 메서드를 사용해서 QDateTime 객체로 변환.
      dtStart = QDateTime.fromMSecsSinceEpoch(int(pvs[0].x()))
      dtLast = dtStart.addSecs(self.last_viewLimit)

      ax = self.EulerChart.axisX()
      ax.setRange(dtStart, dtLast)

      ay = self.EulerChart.axisY()
      ay.setRange(-40, 35)
    except:
      pass

  @pyqtSlot(np.ndarray)
  def Current_GV(self, dotgamma):
    try:
      self.lineEdit_22.setText(str(round(dotgamma[0], 2)))
      self.lineEdit_23.setText(str(round(dotgamma[1], 2)))
      self.lineEdit_24.setText(str(round(dotgamma[2], 2)))
      self.lineEdit_25.setText(str(round(dotgamma[3], 2)))
    except:
      pass

  @pyqtSlot(np.ndarray)
  def Current_GA(self, gamma):
    try:
      self.lineEdit_18.setText(str(round(math.degrees(gamma[0]), 2)))
      self.lineEdit_19.setText(str(round(math.degrees(gamma[1]), 2)))
      self.lineEdit_20.setText(str(round(math.degrees(gamma[2]), 2)))
      self.lineEdit_21.setText(str(round(math.degrees(gamma[3]), 2)))
    except:
      pass

  def btn_clicked_2(self):
    time.sleep(5)
    print('자세제어 시작')
    self.cmg.start()
    self.cmg.finished.connect(self.Current_RPY)
    self.cmg.finished2.connect(self.Current_GV)
    self.cmg.finished3.connect(self.Current_GA)

app = QApplication(sys.argv)
window = MyWindow()
window.show()
app.exec_()