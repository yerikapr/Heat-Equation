class Heat1DSolvers
{
    #length;
    #totalData;
    #partX;
    #diffAr;
    
    constructor(domain, boundary, deltaX, deltaT, diffusivity, initial)
    {
        this.domain = domain;
        this.deltaX = deltaX;
        this.deltaT = deltaT;
        this.diffusivity = diffusivity;
        this.initial = initial;
        this.boundary = boundary;
        this.#length = domain[1] - domain[0];
        this.#totalData = this.#length/this.deltaX + 1;
        this.partX = Array.from({length : this.#totalData}, (_, i) => this.domain[0] + i*deltaX);
        this.#diffAr = this.partX.map(this.diffusivity);

        // console.log('ini this.partX', this.partX);
        // console.log('ini domain', this.domain);
        // console.log('ini deltaX', this.deltaX);
        // console.log('ini initial', this.initial);
    }

    Initializer()
    {
        return this.partX.map(this.initial);
    }

    UpdateDiffusivity()
    {
        this.#diffAr = this.partX.map(this.diffusivity);
    }

    DirichletSolver(lastTemp)
    {
        if(lastTemp == undefined || lastTemp.length == 0)
        {
            return  this.partX.map(this.initial);
        }

        // console.log(this.#diffAr);
        let kmax = Math.max(...this.#diffAr);
        let alpha = this.deltaT*kmax/(this.deltaX**2);
        let temp = [];
        temp[0] = this.boundary[0];
        temp[this.#totalData-1] = this.boundary[1];

        for(let i = 1; i<this.#totalData-1; i++)
        {
           
            temp[i] = lastTemp[i] + alpha*((this.#diffAr[i+1] + this.#diffAr[i])*lastTemp[i+1]/(2*kmax) 
                        - (this.#diffAr[i+1] + 2*this.#diffAr[i] + this.#diffAr[i-1])*lastTemp[i]/(2*kmax)
                        + (this.#diffAr[i] + this.#diffAr[i-1])*lastTemp[i-1]/(2*kmax) );
        }

        return temp;

    }

    NeumannSolver(lastTemp)
    {
        if(lastTemp == undefined || lastTemp.length == 0)
        {
            return  this.partX.map(this.initial);
        }

        let kmax = Math.max(...this.#diffAr);
        let alpha = this.deltaT*kmax/(this.deltaX**2);
        let temp = [];

        temp[0] = lastTemp[0] + alpha*((this.#diffAr[1] + 2*this.#diffAr[0])*(lastTemp[1] - lastTemp[0]) /(2*kmax) 
            - 2*this.#diffAr[0]*this.boundary[0]*this.deltaX /(2*kmax)  ) ;

        temp[this.#totalData-1] = lastTemp[this.#totalData-1] + alpha*( 2*this.#diffAr[this.#totalData-1]*this.deltaX*this.boundary[1] /(2*kmax) -  
            (2*this.#diffAr[this.#totalData-1] + this.#diffAr[this.#totalData-2])*(lastTemp[this.#totalData-1] - lastTemp[this.#totalData-2])/(2*kmax) );

        for(let i = 1; i<this.#totalData-1; i++)
        {
            temp[i] = lastTemp[i] + alpha*((this.#diffAr[i+1] + this.#diffAr[i])*lastTemp[i+1]/(2*kmax) 
                        - (this.#diffAr[i+1] + 2*this.#diffAr[i] + this.#diffAr[i-1])*lastTemp[i]/(2*kmax)
                        + (this.#diffAr[i] + this.#diffAr[i-1])*lastTemp[i-1]/(2*kmax) ); 
        }

        return temp;
    }

}