
const byte interruptPin2 = 2;
const byte interruptPin3 = 3;

volatile byte state1 = false;
volatile byte state2 = false;

unsigned long microsIni = 0;
unsigned long microsFin = 0;
unsigned long microsDif = 0;

void setup() {
  pinMode(interruptPin2, INPUT);
  pinMode(interruptPin3, INPUT);

  attachInterrupt(digitalPinToInterrupt(interruptPin2), int2, FALLING);
  attachInterrupt(digitalPinToInterrupt(interruptPin3), int3, FALLING);

  Serial.begin(115200);
  Serial.println("Serial ok!");
}

void loop(){

  if(state1==true){
	microsIni = micros();
	Serial.print("microsIni:");
	Serial.println(microsIni);
	state1 = false;
  }
 
  if(state2==true){
	microsFin = micros();
	Serial.print("microsFin:");
	Serial.println(microsFin);

	microsDif = microsFin - microsIni;
	Serial.print("microsDif:");
	Serial.println(microsDif);
	state2 = false;
  }
 
}

void int2() {
  state1 = true;
}

void int3() {
  state2 = true;
}
