//- JavaScript source code

//- usm.js ~~
//
//  This is an implementation of Jonas Almeida's "Universal Sequence Maps". The
//  motivation here is to improve the performance of the original version and
//  to add an asynchronous API that enables for distributed parallel execution.
//
//                                                      ~~ (c) SRW, 06 Feb 2013
//                                                  ~~ last updated 07 Feb 2013

(function (global) {
    'use strict';

 // Pragmas

    /*jslint indent: 4, maxlen: 80 */

 // Prerequisites

    // (check for jmat)

    // (check for Quanah)

 // Declarations

    var max2, USM, usm;

 // Definitions

    max2 = function (x) {
     // returns maximum value of array and its index, i.e.  [max,i]
        var xx, y;
        if (Array.isArray(x[0])) {
         // coded only up to 2 dimensions
            xx = x.map(function (xi) {
                return max2(xi);
            }).transpose();
            y = max2(xx[0]);
            return [y[0], [y[1], xx[1][y[1]]]];
        }
        return x.map(function (xi, i) {
            return [xi, i];
        }).reduce(function (a, b) {
            return (a[0] > b[0]) ? a : b;
        });
    };

    USM = function USM(seq, abc, pack, seed) {
     // This function needs documentation.
        var that = this;
        if (seq) {
         // find out if this is a sequence or the url of a fastA file with one
            if (seq.length > 15) { // it could be a url
                if (!!seq.slice(0, 10).match(/:\/\//)) {
                    global.console.log('using proxy ' + global.jmat.webrwUrl +
                            '\nto get fastA file ' + seq + ' ...');
                    that.loadFasta(seq, function (x) {
                     // This function needs documentation.
                        global.console.log('... file loaded ...');
                    }, abc);
                } else {
                    that.encode(seq, abc, pack);
                }
            } else {
                that.encode(seq, abc, pack, seed);
            }
        } else {
            if (abc) {
             // if alphabet is provided then identify cube anyway to enable
             // decoding
                that.encode(abc, abc, pack, seed);
            }
        }
        return that;
    };

    usm = function (seq, abc, pack, seed) {
     // Universal Sequence Map
        return new USM(seq, abc, pack, seed);
    };

 // Prototype definitions

    USM.prototype.align = function (sprobe, sbase) {
     // align sequence sprobe to this sequence
        if (typeof sprobe === 'string') {
         // in case the actual sequence was inputed
            sprobe = usm(sprobe, this.abc, this.pack);
        }
        if (!sbase) {
            sbase = this;
        }
     // start by considering a complete match and zoom in until such a subset
     // is found
        var A, n, nn;
        A = {
            posBase: [0],
            posProbe: [0],
            match: [0],
            ind: 0
        };
        n = sprobe.seq.length;
        nn = sbase.seq.length;
        this.alignUsm = function (posStart, posEnd) {
         // defined here to capture align closure
            global.console.log([posStart, posEnd]);
            var d, dF, dFi, i, inc, j, mm, res;
            inc = Math.floor((posEnd - posStart) / 2); // increment
            i = posStart + inc;
            res = 32; //some resolution
            d = sbase.cgrForward.map(function (cf, ii) {
                return sbase.distCGR(cf, sprobe.cgrForward[i]) +
                        sbase.distCGR(sbase.cgrBackward[ii],
                                sprobe.cgrBackward[i]);
            });
            mm = max2(d); // I removed the jmat dependency here by inlining
            if (mm[0] > 0) {
                mm[0] -= 1;
            }
            j = A.posProbe.length;
         // forward distance
            dF = sbase.distCGR(sbase.cgrForward[mm[1]], sprobe.cgrForward[i]);
            if (dF > 0) {
                dF -= 1;
            }
         // using foward CGR to find start position
            A.posBase[j] = mm[1] - dF;
            if (A.posBase[j] < 0) {
             // check lower boundary for Base
                mm[0] -= (dF - mm[1]);
                dF = mm[1];
                A.posBase[j] = 0;
            }
            if ((A.posBase[j] + dF) > nn) {
             // check upper boundary for Base
                mm[0] -= A.posBase[j] + mm[0] - nn;
            }
            A.posProbe[j] = i - dF;
            if (A.posProbe[j] < 0) {
             // check lower boundary for Probe
                mm[0] -= (dF - i);
                dF = i;
                A.posProbe[j] = 0;
                A.posBase[j] = mm[1] - dF;
            }
            if ((A.posProbe[j] + mm[0]) > n) {
             // check upper boundary of Probe
                mm[0] -= A.posProbe[j] + mm[0] - n;
            }
            if (mm[0] > res) {
             // check if numerical resolution might have been exceeded
             // additional forward distance
                dFi = sbase.distCGR(sbase.cgrForward[A.posBase[j]],
                        sprobe.cgrForward[A.posProbe[j]]) - 1;
                while (dFi > 0) { // if some was found
                    A.posProbe[j] -= dFi;
                    A.posBase[j] -= dFi;
                    mm[0] += dFi;
                    dFi = sbase.distCGR(sbase.cgrForward[A.posBase[j]],
                            sprobe.cgrForward[A.posProbe[j]]) - 1;
                }
             // backward distance at the end of match segment
                dFi = sbase.distCGR(sbase.cgrBackward[A.posBase[j] +
                        mm[0] - 1], sprobe.cgrBackward[A.posProbe[j] +
                        mm[0] - 1]) - 1;
                while (dFi > 0) {
                 // if some was found
                    mm[0] += dFi;
                    dFi = sbase.distCGR(sbase.cgrBackward[A.posBase[j] +
                            mm[0] - 1], sprobe.cgrBackward[A.posProbe[j] +
                            mm[0] - 1]) - 1;
                }
            }
            A.match[j] = mm[0];
            if (mm[0] > (A.match[A.ind])) {
             // if this is the best match yet
                A.ind = j;
            }
            if (A.match[A.ind] < inc) {
             // if best match is shorter than increment, zoom into its halves
                inc = Math.floor(inc / 2);
                this.alignUsm(posStart, inc);
                this.alignUsm(inc + 1, posEnd);
            }
        };
        this.alignUsm(0, n);
        global.console.log('largest identical segment has length ' +
                A.match[A.ind] + ' and aligns with position ' +
                A.posBase[A.ind] + ' in base sequence and position ' +
                A.posProbe[A.ind] + ' in probe sequence');
        return A;
    };

    USM.prototype.alignQ = function (sprobe, sbase) {
     // This function provides an asynchronous wrapper for the usual
     // `this.align` method using Quanah (http://wilkinson.github.com/quanah/).
        if (Object.prototype.hasOwnProperty('Q') === false) {
            throw new Error('Quanah is not loaded.');
        }
        var Q, that, y;
        Q = Object.prototype.Q;
        that = this;                    //- the current USM object
        y = Q.avar();
        y.Q(function (evt) {
         // This function needs documentation.
            if (that.hasOwnProperty('cgrBackward') === false) {
                global.setTimeout(y.revive, 0);
                return evt.stay('Waiting for indexing to finish ...');
            }
            y.val = that.align(sprobe, sbase);
            return evt.exit();
        });
        return y;
    };

    USM.prototype.alphabet = function (seqString) {
     // extracts alphabet
        var i, that;
        that = this;
        if (!seqString) {
         // uses own string if argument not provided
            seqString = that.seq;
        }
        that.abc = '';
        for (i = 0; i < seqString.length; i += 1) {
            if (!that.abc.match(new RegExp(seqString[i]))) {
                that.abc += seqString[i];
            }
        }
        return that.abc.sort(); // using overloaded String.sort()
    };

    USM.prototype.bin2int = function (B) {
     // converts binary vector into integer
        return B.slice().reverse().map(function (x, i) {
            return x * Math.pow(2, i);
        }).reduce(function (a, b) {
            return a + b;
        });
    };

    USM.prototype.cgr = function (bin, s, ith) {
     // CGR with recursive seed
        if (!ith) {
            ith = 0;
        }
        var i, n, y;
        n = bin.length;
        if (!s) {
         // start seed with last value of bin
            s = bin[bin.length - 1];
        }
        y = [];
        y[0] = s + ((bin[0] - s) / 2);
        for (i = 1; i < n; i += 1) {
            y[i] = y[i - 1] + ((bin[i] - y[i - 1]) / 2);
        }
     // check recursive seed
        if ((s !== y[y.length - 1]) && (ith < 64)) {
            y = this.cgr(bin, y[y.length - 1], ith + 1);
        }
        return y;
    };

    USM.prototype.cgr2 = function (bin, seed) {
     // CGR with 1/2 seed
        var i, y;
        y = [(1 / 2) + ((bin[0] - (1 / 2)) / 2)];
        for (i = 1; i < bin.length; i += 1) {
            y[i] = y[i - 1] + ((bin[i] - y[i - 1]) / 2);
        }
        return y;
    };

    USM.prototype.cgrLong = function (ii, direction, s) {
     // one dimension at a time
        var bin, c, i, n, that;
        that = this;
        bin = that.bin[ii];
        n = bin.length; // binning threads
        if (direction === 'cgrBackward') {
            bin.reverse();
            c = that.cgrBackward[ii];
        } else {
            c = that.cgrForward[ii];
        }
        if (!s) {
         // seed
            s = Math.random();
            for (i = (n - 128); i < n; i += 1) {
                s = s + (bin[i] - s) / 2;
            }
        }
        c[0] = s + ((bin[0] - s) / 2);
        for (i = 1; i < n; i += 1) {
            c[i] = c[i - 1] + ((bin[i] - c[i - 1]) / 2);
        }
        return 'done';
    };

    USM.prototype.decode = function (xy, n) {
     // decode numerical coordinates
        var abc, bb, bin2int, decodeBin;
        abc = this.abc;
        bin2int = this.bin2int;
        decodeBin = this.decodeBin;
        bb = xy.map(decodeBin);
        bb = this.transpose(bb).map(bin2int);
        bb = bb.map(function (x) {
            return abc[x];
        });
        return bb.toString().replace(/,/g, '');
    };

    USM.prototype.decodeBin = function (x, n) {
     // decompose single coordinate, x, into a binary array of length <= n
        if (!n) {
            n = Infinity;
        }
        var i, y, z;
        i = 0;
        y = [];
        z = [];
        x *= 2;
     // just in case x is 0, 1 or even 1/2 (impossible) y is still populated
        y[0] = x;
        while ((x !== Math.round(x)) && (i < n)) {
            if (x > 1) {
                y[i] = 1;
                x -= 1;
            } else {
                y[i] = 0;
            }
            x *= 2;
            i += 1;
        }
        return y;
    };

    USM.prototype.dist = function (x, y) {
     // Equation 5: distance between two usm coordinates for a position, i
     // each provided as a two element Array [cgrForward[i],cgrBackward[i]]
        var d = this.distCGR(x[0], y[0]) + this.distCGR(x[1], y[1]) - 1;
        return (d < 0) ? 0 : d;
    };

    USM.prototype.distCGR = function (a, b) {
     // distance between two CGR positions, a and b
        var dist = this.L;
        return this.transpose([a, b]).map(function (x) {
            return dist(x[0], x[1]);
        }).min();
    };

    USM.prototype.distMap = function (sprobe, sbase) {
     // Distance Map to a new probing sequence, sprobe
        if (typeof sprobe === 'string') {
         // in case the actual sequence was specified
            sprobe = usm(sprobe, this.abc, this.pack);
        }
        if (!sbase) {
            sbase = this;
        }
        return sbase.usm.map(function (x) {
            return sprobe.usm.map(function (y) {
                return sbase.dist(x, y);
            });
        });
    };

    USM.prototype.distProfile = function (s) {
     // Distance to another sequence, s
        if (typeof s === 'string') {
            s = usm(s, this.abc, this.pack);
        }
        if (s.abc !== this.abc) {
            throw new Error('unequal alphabets');
        }
        if (s.pack !== this.pack) {
            throw new Error('unequal encode packing');
        }
        var b, f;
        f = this.cgrForward.map(function (x) {
            var x0 = x;
            return s.cgrForward.map(function (x) {
                return s.distCGR(x0, x);
            }).sum();
        });
        b = this.cgrBackward.map(function (x) {
            var x0 = x;
            return s.cgrBackward.map(function (x) {
                return s.distCGR(x0, x);
            }).sum();
        });
        return [f, b].transpose().sum();
    };

    USM.prototype.encode = function (seq, abc, pack, seed) {
        if (!seq) {
            seq = this.seq;
        } else {
            this.seq = seq;
        }
        if (!this.seq) {
            throw new Error('Sequence not provided');
        }
        if (abc) {
            if (abc.length === 0) {
             // in case `abc` === `''`
                abc = undefined;
            }
        }
        if (abc) {
            this.abc = abc;
        }
        if (!this.abc) {
            this.abc = this.alphabet();
        }
        var i, m, n;
        m = this.abc.length;
        n = this.seq.length;
        this.cube = [];
        this.str2cube(pack);
        this.cgrForward = [];
        this.cgrBackward = [];
        if (!seed) {
            for (i = 0; i < this.bin.length; i += 1) {
                this.cgrForward[i] = this.cgr(this.bin[i]);
                this.cgrBackward[i] = this.cgr(this.bin[i]
                        .slice().reverse()).slice().reverse();
            }
        } else {
         // seeded CGR
            for (i = 0; i < this.bin.length; i += 1) {
                this.cgrForward[i] = this.cgr2(this.bin[i], seed);
                this.cgrBackward[i] = this.cgr2(this.bin[i]
                        .slice().reverse(), seed).slice().reverse();
            }
        }
        this.cgrForward = this.transpose(this.cgrForward);
        this.cgrBackward = this.transpose(this.cgrBackward);
        this.usm = [];
        for (i = 0; i < n; i += 1) {
         // In a serious application .usm is all we'd need to keep .cgr--- etc
         // could all be deleted
            this.usm[i] = [this.cgrForward[i], this.cgrBackward[i]];
        }
        return;
    };

    USM.prototype.encodeLong = function (seq, abc, pack, seed) {
     // encoding long sequences by writting directly to the usm instance
        var c, i, j, m, n, that;
        c = [];
        that = this;
        if (!!seq) {
            that.seq = seq;
        }
        if (!!abc) {
            that.abc = abc;
        }
        if (!that.seq) {
            throw new Error('Sequence not provided');
        }
        if (!that.abc) {
            global.console.log('find alphabet ...');
            that.abc = that.alphabet();
        }
        global.console.log('... alphabet: ' + that.abc);
        that.cube = [];
        global.console.log('packing USM space ...');
        that.str2cube(pack, true);
        global.console.log('...', that.cube);
        global.console.log('filling USM space ...');
        m = that.cube.length;
        n = that.seq.length;
        that.cgrForward = [];
        that.cgrBackward = [];
        for (i = 0; i < m; i += 1) {
            that.cgrForward[i] = [];
            that.cgrBackward[i] = [];
        }
        for (i = 0; i < m; i += 1) {
            global.console.log(((i * 2) + 1) + '/' + (that.cube.length * 2) +
                    ' mapping axis <' + that.cube[i] + '> forward');
            that.cgrLong(i, 'cgrForward');
            global.console.log(((i * 2) + 2) + '/' + (that.cube.length * 2) +
                    ' mapping axis <' + that.cube[i] + '> backward');
            that.cgrLong(i, 'cgrBackward');
            that.bin[i] = []; // free memory
        }
        global.console.log('packing CGR coordinates ...');
        global.console.log('... forward ...');
        for (i = 0; i < n; i += 1) {
            c[i] = [];
            for (j = 0; j < m; j += 1) {
                c[i][j] = this.cgrForward[j][i];
            }
        }
        that.cgrForward = c;
        global.console.log('... backward ...');
        c = [];
        for (i = 0; i < n; i += 1) {
            c[i] = [];
            for (j = 0; j < m; j += 1) {
                c[i][j] = that.cgrBackward[j][i];
            }
        }
        that.cgrBackward = c.reverse();
        c = [];
        global.console.log('USMapping done');
        return;
    };

    USM.prototype.L = function (a, b) {
     // distance between two coordinates
        var d, temp;
        d = 0;
        temp = 1;
        while ((temp !== Infinity) &&
                (Math.round(a * temp) === Math.round(b * temp))) {
            d += 1;
            temp *= 2;
        }
        if (temp === Infinity) {
         // stop at the numerical resolution
            d = 64;
        }
        return d;
    };

    USM.prototype.loadFasta = function (url, callback, abc) {
     // load long sequences from a FastA file
        var that = this;
        if (!!abc) {
            that.abc = abc;
        }
        global.jmat.get(url, function (x) {
         // encode sequences
            if (!!callback) {
             // in case this function was called with a callback do that first
                callback(x);
            }
            var n, nn;
         // let's encode it now
         //
         // fastA head identification, everything else should be a long sequence
            global.console.log(x[0]);
            n = x.length;
            that.seq = x.splice(1, n).join('');
            x = ''; // to free memory
            nn = that.seq.length;
            global.console.log('... found ' + nn + ' units divided in ' + n +
                    ' segments,');
            that.encodeLong(); // the seq will be picked from `that.seq`
            return;
        });
        return 'using proxy ' + global.jmat.webrwUrl + ' to get fastA file ' +
                url + ' ...';
    };

    USM.prototype.mapReduce = function (x, map, reduce) {
     // This function needs documentation.
        return reduce(x.map(map));
    };

    USM.prototype.str2cube = function (pack, show) {
        var abc, f, i, j, L, m, mm, n;
        m = this.abc.length;
        n = this.seq.length;
        if (!pack) {
         // default packing method
            pack = 'compact';
        }
        this.pack = pack;
        this.bin = [];
        switch (pack) {
        case 'sparse':
            f = function (si) {
             // This function needs documentation.
                return (si === this.abc[j]) ? 0 : 1;
            };
            for (j = 0; j < m; j += 1) {
                this.cube[j] = this.abc[j];
                this.bin[j] = this.seq.split('').map(f);
            }
            break;
        case 'compact':
            L = Math.ceil(Math.log(m) / Math.log(2)); // map dimension
            mm = Math.pow(2, L); // maximum length of this alphabet
            f = function (si) {
             // This function needs documentation.
                return (abc.match(new RegExp(si))) ? 0 : 1;
            };
            for (j = 0; j < L; j += 1) {
                abc = '';
                mm = mm / 2;
                for (i = 0; i < m; i += (mm * 2)) {
                    abc += this.abc.slice(i, i + mm);
                }
                this.cube[j] = abc;
                if (show) {
                    global.console.log((j + 1) + '/' + L + ' filling axis <' +
                            abc + '>');
                }
                this.bin[j] = this.seq.split('').map(f);
            }
            break;
        default:
         // (placeholder)
        }
        return;
    };

    USM.prototype.transpose = function (m) {
        var i, j, t;
        t = [];
        for (i = 0; i < m[0].length; i += 1) {
            t[i] = [];
            for (j = 0; j < m.length; j += 1) {
                t[i][j] = m[j][i];
            }
        }
        return t;
    };

    USM.prototype.transpose2 = function (M) {
     // I'm honestly not sure if this function is ever used ...
        M[0].map(function (mi, i) {
            var MM = [];
        });
        return;
    };

 // Out-of-scope definitions

    Array.prototype.max = function () {
     // returns maximum value
        return this.reduce(function (x1, x2) {
         // This function needs documentation.
            return (x1 > x2) ? x1 : x2;
        });
    };

    Array.prototype.min = function () {
     // returns maximum value
        return this.reduce(function (x1, x2) {
         // This function needs documentation.
            return (x1 < x2) ? x1 : x2;
        });
    };

    Array.prototype.sum = function () {
     // returns sum of all values
        return this.reduce(function (x1, x2) {
         // This function needs documentation.
            return x1 + x2;
        });
    };

    Array.prototype.transpose = function () {
     // written for arrays of arrays
        var i, j, m, t;
        if (!Array.isArray(this[0])) {
            m = this.map(function (x) {
             // This function needs documentation.
                return [x];
            });
        } else {
            m = this;
        }
        t = [];
        for (i = 0; i < m[0].length; i += 1) {
            t[i] = [];
            for (j = 0; j < m.length; j += 1) {
                t[i][j] = m[j][i];
            }
        }
        if ((Array.isArray(t[0])) && (t[0].length === 1)) {
            t = t.map(function (x) {
             // This function needs documentation.
                return x[0];
            });
        }
        return t;
    };

    String.prototype.reverse = function () {
     // sort characters of a string
        return this.split('').reverse().toString().replace(/,/g, '');
    };

    String.prototype.sort = function () {
     // sort characters of a string
        return this.split('').sort().toString().replace(/,/g, '');
    };

    global.usm = usm;

 // Invocations

    (function () {
     // This function needs documentation.
        /*jslint browser: true */
        if (global.hasOwnProperty('document') === false) {
            return;
        }
        if (global.hasOwnProperty('jmat')) {
            return;
        }
        var s = global.document.createElement('script');
        s.src = 'https://jmat.googlecode.com/git/jmat.js';
        global.document.head.appendChild(s);
        return;
    }());

    (function () {
     // This function needs documentation.
        /*jslint node: true */
        if (typeof module === 'object') {
            module.exports = usm;
        }
        return;
    }());

    global.console.log('CGR toolbox :-)');

 // That's all, folks!

    return;

}(Function.prototype.call.call(function (that) {
    'use strict';

 // This strict anonymous closure encapsulates the logic for detecting which
 // object in the environment should be treated as _the_ global object. It's
 // not as easy as you may think -- strict mode disables the `call` method's
 // default behavior of replacing `null` with the global object. Luckily, we
 // can work around that by passing a reference to the enclosing scope as an
 // argument at the same time and testing to see if strict mode has done its
 // deed. This task is not hard in the usual browser context because we know
 // that the global object is `window`, but CommonJS implementations such as
 // RingoJS confound the issue by modifying the scope chain, running scripts
 // in sandboxed contexts, and using identifiers like `global` carelessly ...

    /*jslint indent: 4, maxlen: 80 */
    /*global global: true */
    /*properties global */

    if (this === null) {

     // Strict mode has captured us, but we already passed a reference :-)

        return (typeof global === 'object') ? global : that;

    }

 // Strict mode isn't supported in this environment, but we need to make sure
 // we don't get fooled by Rhino's `global` function.

    return (typeof this.global === 'object') ? this.global : this;

}, null, this)));

//- vim:set syntax=javascript:
