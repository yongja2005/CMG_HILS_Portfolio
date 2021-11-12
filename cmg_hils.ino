int PIN_dir = 12;
int PIN_step = 13;
const int encoderPinA = 2;
const int encoderPinB = 3;
String incommingByte = "";
int encoderPos = 0;
#define radians(deg) ((deg)*DEG_TO_RAD)
const float ratio = radians(360./4000.);

void doEncoderA(){  
  encoderPos += (digitalRead(encoderPinA)==digitalRead(encoderPinB))?1:-1;
}
void doEncoderB(){  
  encoderPos += (digitalRead(encoderPinA)==digitalRead(encoderPinB))?-1:1;
}

void setup(){
  pinMode(encoderPinA, INPUT_PULLUP);
  attachInterrupt(0, doEncoderA, CHANGE);
 
  pinMode(encoderPinB, INPUT_PULLUP);
  attachInterrupt(1, doEncoderB, CHANGE);
  
  Serial.begin(115200);
  
  pinMode(PIN_dir,OUTPUT);
  pinMode(PIN_step,OUTPUT);
  pinMode(11, OUTPUT);
  digitalWrite(11, HIGH);
  pinMode(6, OUTPUT);
  pinMode(7, OUTPUT);
}
 
void loop(){ 
  
  if (Serial.available() > 0){
  digitalWrite(6, HIGH); 
  incommingByte = Serial.readStringUntil('Q');
  int index1 = incommingByte.indexOf(',');
  int index2 = incommingByte.length();
  if (incommingByte.substring(0, index1).indexOf('-') == 0) {
    int s = incommingByte.substring(0, index1).toInt() * -1;
  }
  int s = incommingByte.substring(0, index1).toInt();
  int ss = incommingByte.substring(index1+1, index2).toInt();
  int step_size = abs(s)/0.05625;
  int step_delay = 0.1/step_size*1000000;
  int sss = 0.1/s*1000000;
   for(int i = 0; i < abs(s); i++){// 1회 동작 당 값
      if (s > 0) {
        digitalWrite(PIN_dir,HIGH);
      }
      else {
        digitalWrite(PIN_dir,LOW);
      }
      digitalWrite(PIN_step,HIGH);    
      delayMicroseconds(abs(sss/2));    //회전 속도
      digitalWrite(PIN_step,LOW);
      delayMicroseconds(abs(sss/2));    
    }
    Serial.println(float(encoderPos)*ratio);
  }
}
