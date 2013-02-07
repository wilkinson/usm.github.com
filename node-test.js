//- JavaScript source code

//- node-test.js ~~
//
//  The purpose of this file is NOT to test USM for use as in Node.js, because
//  server-side programming is going extinct! This just helps me test that the
//  new implementation of USM does not depend on any functionality outside what
//  is available in the language itself.
//
//                                                      ~~ (c) SRW, 06 Feb 2013
//                                                  ~~ last updated 07 Feb 2013

(function () {
    'use strict';

 // Pragmas

    /*jslint indent: 4, maxlen: 80, node: true */

 // Declarations

    var usm;

 // Definitions

    usm = require('./usm');

 // Demonstrations

    (function () {
        /*jslint newcap: true */
        console.log(new usm('CAT').align('GAT'));
        return;
    }());

 // That's all, folks!

    return;

}());

//- vim:set syntax=javascript:
