const xSliderTolerance = 5;
const diffusivityRange = 10;
const Tmin = 0;
const Tmax = 100;

class Menu1 {
    constructor() {
        this.temp = []; //temperature
        this.windowWidth = 500;
        this.windowHeight = 300;
        this.windowSize = [this.windowWidth, this.windowHeight];
        this.partition = 100;
        this.deltaX = this.windowSize[0] / this.partition;
        this.deltaT = 5;
        this.windowPos = [150, height / 2 - this.windowHeight / 2 - 100];
        this.windowTopLeft = [this.windowPos[0], this.windowPos[1]];
        this.windowTopRight = [this.windowPos[0] + this.windowWidth, this.windowPos[1]];
        this.windowBottomLeft = [this.windowPos[0], this.windowPos[1] + this.windowHeight];
        this.windowBottomRight = [this.windowPos[0] + this.windowWidth, this.windowPos[1] + this.windowHeight];

        this.formulaPanel = new FormulaPanel(this);
        this.temperaturePanel = new TemperaturePanel({topY: this.windowTopLeft[1]});
        this.timePanel = new TimePanel(this);
        this.distributionPanel = new DistributionPanel(this);
        this.graphicPanel = new GraphicPanel(this);
        
        this.isPlaying = false;
        this.isRunning = false;
        
        this.DiffusivityPanelA();
        this.DiffusivityPanelB();
        this.ButtonPanel();
        this.formulaPanel.create();

        this.formula = new Heat1DSolvers([this.windowPos[0], this.windowPos[0] + this.windowWidth], [0,100], this.deltaX, this.deltaT, this.diffusivity.bind(this), this.initialDistribution.bind(this));
    }
 
    display() {
        this.WindowPanel();
        this.timePanel.display(this);
        this.temperaturePanel.display();
        this.diffusivityPanelA.display();
        this.diffusivityPanelB.display();
        this.distributionPanel.display();
        this.graphicPanel.display(this);
        
        if (this.isRunning) {
            this.distributionPanel.disabled();
            this.diffusivityPanelA.disabled();
            this.diffusivityPanelB.disabled();
            this.formulaPanel.disable();
            this.formula.UpdateDiffusivity();

            this.graphicPanel.drawGraph(this);
        } else {
            this.distributionPanel.enabled();
            this.diffusivityPanelA.enabled();
            this.diffusivityPanelB.enabled();
            this.temp = this.formula.Initializer();
            this.graphicPanel.drawGraph(this);
        }

        if (this.isPlaying) {
            this.playButton.attribute('disabled', '');
            this.playButton.style('background-color', '#cccccc');
            this.playButton.style('color', '#666666');

            this.pauseButton.removeAttribute('disabled');
            this.pauseButton.style('background-color', '#ffc107');
            this.pauseButton.style('color', '#ffffff');

            if (this.formulaPanel.selected === 1) {
                this.formula.boundary = [0,100];
                this.temp = this.formula.DirichletSolver(this.temp);
            } else {
                this.formula.boundary = [0,0];
                this.temp = this.formula.NeumannSolver(this.temp);
            }


            // console.log(`cek DIFUSITI A`, this.diffusivityPanelA.getValue());
            // console.log(`cek DIFUSITI B`, this.diffusivityPanelB.getValue());
        } else {
            this.playButton.removeAttribute('disabled');
            this.playButton.style('background-color', '#4CAF50');
            this.playButton.style('color', '#ffffff');

            this.pauseButton.attribute('disabled', '');
            this.pauseButton.style('background-color', '#cccccc');
            this.pauseButton.style('color', '#666666');
        }
    }

    WindowPanel() {
        colorMode(RGB, 255);
        noStroke();
        for (let i = 0; i < this.windowSize[0] / this.deltaX; i++) {
            fill(this.tempToColor(this.temp[i]));
            rect(this.windowPos[0] + i * this.deltaX , this.windowPos[1], this.deltaX, this.windowSize[1]);
        }
        stroke(0, 0, 0);
        strokeWeight(1);
    }

    DiffusivityPanelA() {
        let posY = this.distributionPanel.getBottomLeft()[1];
        this.diffusivityPanelA = new Diffusivity({
            props: this,
            label: 'A',
            x: this.windowBottomLeft[0],
            y: posY,
        }, 10);
        this.diffusivityPanelA.create();
    }

    DiffusivityPanelB() {
        let posX = this.windowBottomLeft[0] + (this.windowWidth / 2) + 15;
        let posY = this.distributionPanel.getBottomLeft()[1];
        this.diffusivityPanelB = new Diffusivity({
            props: this,
            label: 'B',
            x: posX,
            y: posY,
        });
        this.diffusivityPanelB.create();
    }

    ButtonPanel() {
        let posY = this.diffusivityPanelA.getBottomLeft()[1] + 60;
        let posX = this.windowBottomLeft[0];

        this.playButton = createButton("Play");
        this.playButton.position(posX, posY);
        this.playButton.style('background-color', '#4CAF50');
        this.playButton.style('width', '90px');
        this.playButton.mouseClicked(() => {
            this.timePanel.start();
            this.isPlaying = true;
            this.isRunning = true;
        });

        posX = posX + 100
        this.pauseButton = createButton("Pause");
        this.pauseButton.position(posX, posY);
        this.pauseButton.style('width', '90px');
        this.pauseButton.style('background-color', '#ffc107');
        this.pauseButton.mouseClicked(() => {
            this.isPlaying = false;
            this.timePanel.pause();
        });

        posX = posX + 100
        this.resetButton = createButton("Reset");
        this.resetButton.position(posX, posY);
        this.resetButton.style('width', '90px');
        this.resetButton.style('background-color', '#dc3545');
        this.resetButton.mousePressed(() => {
            // this.temp = [];
            // redraw(1);
            window.location.reload();
        });
    }

    tempToColor(temperature) {
        let deltaT = Tmax - Tmin;
        if (temperature == null) return [255, 255, 255];
        if (deltaT === 0) return [255, 255, 255]; // Menghindari pembagian dengan nol
        let r = 255 * (temperature - Tmin) / deltaT;
        let g = 0;
        let b = 255 * (Tmax - temperature) / deltaT;
        return [r, g , b];
    }

    diffusivity(x) {
        let posX = x - this.windowPos[0];
        let val1 = this.diffusivityPanelA.getValue();
        let val2 = this.diffusivityPanelB.getValue();
        //let val1 = 1.0;
        //let val2 = 10;

        console.log('====== CEK A', val1);
        console.log('====== CEK B', val2);
    
        return posX < this.windowWidth ? val1 : val2;
    }

    initialDistribution(x) {
        let posX = x - this.windowPos[0];
        return posX < this.distributionPanel.getValue() ? 0 : 100;
        // return posX < this.windowWidth / 2 ? 0 : 100;
        // return 100 * (x - this.windowPos[0])/this.windowWidth;
    }
}

class FormulaPanel {
    constructor(props) {
        this.dropdown = null;
        this.windowTopRight = props.windowTopRight;
        this.selected = 1;
    }

    create() {
        this.dropdown = createSelect();
        this.dropdown.position(this.windowTopRight[0] - 105, this.windowTopRight[1] - 20);
        this.dropdown.option('Dirichlet Solver');
        this.dropdown.option('Neumann Solver');
        this.dropdown.changed(() => this.onChange());
    }

    onChange() {
        let selectedOption = this.dropdown.value();

        if (selectedOption === 'Dirichlet Solver') {
            this.selected = 1;
        } else if (selectedOption === 'Neumann Solver') {
            this.selected = 2;
        }
    }

    getValue() {
        return this.selected;
    }

    disable() {
        this.dropdown.attribute('disabled', '');
    }

    enable() {
        this.dropdown.removeAttribute('disabled');
    }
}

class TimePanel {
    constructor(props) {
        this.leftX = props.windowTopLeft[0];
        this.leftY = props.windowTopLeft[1] - 20;
        this.startTime = 0;
        this.elapsedTime = 0;
        this.isPlaying = props.isPlaying;
    }

    display(props) {
        this.isPlaying = props.isPlaying;

        // Jika timer sedang berjalan, hitung waktu yang berlalu
        if (this.isPlaying) {
            this.elapsedTime = millis() - this.startTime;
        }
    
         // Hitung jam, menit, detik, dan milidetik
        let hours = Math.floor(this.elapsedTime / (1000 * 60 * 60));
        let minutes = Math.floor((this.elapsedTime % (1000 * 60 * 60)) / (1000 * 60));
        let seconds = Math.floor((this.elapsedTime % (1000 * 60)) / 1000);
        let milliseconds = this.elapsedTime % 1000;

        // Format menjadi dua digit jika dibawah 10
        hours = nf(hours, 2);
        minutes = nf(minutes, 2);
        seconds = nf(seconds, 2);
        milliseconds = nf(milliseconds, 3, 0).replace(".", "");

        fill(0);
        textSize(20);
        text(`Waktu: ${hours}:${minutes}:${seconds}:${milliseconds}`, this.leftX, this.leftY);
        textSize(12);
    }

    start() {
        if (!this.isPlaying) {
          this.startTime = millis() - this.elapsedTime; // Mulai dari waktu yang tersimpan
        }
      }
    
    pause() {
        if (this.isPlaying) {
            this.elapsedTime = millis() - this.startTime; // Simpan waktu yang berlalu
        }                   
    }
}

class TemperaturePanel {
    constructor({topY = topY}) {
        this.panelWidth = 20;
        this.panelHeight = 300;
        this.topY = topY;
        this.bottomY = this.topY + this.panelHeight;
        this.leftX = 50;
        this.rightX = this.leftX + this.panelWidth;
        this.partition = 100;
        this.partHeight = this.panelHeight / this.partition;
    }

    display() {
        for (let i = 0; i < this.partition; i++) {
            let y = (i * this.partHeight) + this.topY;
            // Interpolasi warna
            let t = i / this.partition; // Proporsi posisi bagian dari atas ke bawah
            let r = lerp(255, 0, t);
            let g = 0;
            let b = lerp(0, 255, t);

            fill(r, g, b);
            rect(this.leftX, y, this.panelWidth, this.partHeight);
            if ((i) % 20 == 0) {
                line(this.leftX, y, this.rightX + 5, y);
                textAlign(LEFT, CENTER);
                text(this.partition - i, this.rightX + 8, y);
            }
        }

        line(this.leftX, this.bottomY, this.rightX + 5, this.bottomY);
        textAlign(LEFT, CENTER);
        text('0', this.rightX + 8, this.bottomY);

        fill(0,0,0);
        text('Pita Suhu T (ºC)', this.leftX - 25, this.topY - 20);
    }
}