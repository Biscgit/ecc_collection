use std::{arch::asm, cmp::PartialEq, ops};

use rand::Rng;
use integer_sqrt::IntegerSquareRoot;

const MAX_FACTOR: i128 = 1_000;
const MAX_ITERATIONS: u32 = 10_000;


// helper functions
/// returns the most significant bit of a number
#[cfg(target_arch = "x86_64")]
pub fn get_msb_position(number: i128) -> u8 {
    let (high, low): (u64, u64) = ((number.abs() >> 64) as u64, number.abs() as u64);
    let (msb_high, msb_low): (u32, u32);

    unsafe {
        asm!(
        "bsr {result:r}, {input:r}",
        result = lateout(reg) msb_high,
        input = in(reg) high,
        );
    }

    unsafe {
        asm!(
        "bsr {result:r}, {input:r}",
        result = lateout(reg) msb_low,
        input = in(reg) low,
        );
    }

    return match msb_high > 0 {
        true => { msb_high + 32 }
        false => { msb_low }
    } as u8;
}

/// runs the square_and_multiply algorithm for exponentiation
pub fn mod_pow(base: i128, exponent: i128, modulo: i128) -> i128 {
    // get the position of the most significant bit and run algorithm
    let msb = get_msb_position(exponent);
    let mut result: i128 = 1;

    for index in (0..=msb).rev() {
        // square
        result = (result * result) % modulo;

        // multiply
        if (exponent >> index) & 0b1 == 1 {
            result = (result * base) % modulo;
        }
    }

    return result;
}

/// runs the double_and_add algorithm to multiply two numbers
pub fn mod_mul(base: i128, factor: i128, modulo: i128) -> i128 {
    // switch values to increase performance with large factors and small bases
    let (base, factor) = match factor.abs() > base.abs() {
        true => (factor, base),
        false => (base, factor)
    };

    // get the position of the most significant bit and run algorithm
    let msb = get_msb_position(factor);
    let mut result: i128 = 0;

    for index in (0..=msb).rev() {
        // double
        result = (result + result).rem_euclid(modulo);

        // add
        if (factor >> index) & 0b1 == 1 {
            result = (result + base).rem_euclid(modulo);
        }
    }

    return result;
}

/// returns the modular inverse of a number if it exists
pub fn mod_inv(number: i128, modulo: i128) -> Option<i128> {
    let (g, result, _) = euclid_gcd(number.rem_euclid(modulo), modulo);

    match g {
        1 => Some(result.rem_euclid(modulo)),
        _ => None,
    }
}

/// interface for euclidean gcd
pub fn gcd(number1: i128, number2: i128) -> i128 {
    euclid_gcd(number1, number2).0
}

/// executes euclidean gcd
fn euclid_gcd(number1: i128, number2: i128) -> (i128, i128, i128) {
    match number1 {
        0 => { (number2, 0, 1) }
        _ => {
            let (g, x, y) = euclid_gcd(number2.rem_euclid(number1), number1);
            (g, y - (number2 / number1) * x, x)
        }
    }
}


// Lenstra and EC
#[derive(Copy, Clone)]
pub struct WeierStrass {
    a: i128,
    b: i128,
    p: i128,
}

impl WeierStrass {
    pub fn new(a: i128, b: i128, p: i128) -> Option<Self> {
        match (4 * mod_pow(a, 3, p) + 27 * mod_pow(b, 2, p)) % p {
            0 => None,
            _ => Some(WeierStrass { a, b, p })
        }
    }
}

impl PartialEq for WeierStrass {
    fn eq(&self, other: &Self) -> bool {
        self.a == other.a && self.b == other.b && self.p == other.p
    }
}

#[derive(Copy, Clone)]
pub struct WeierStrassPoint {
    x: i128,
    y: i128,
    y_infinite: bool,
    curve: WeierStrass,
}

impl WeierStrassPoint {
    pub fn new(x: i128, y: i128, curve: WeierStrass) -> Self {
        WeierStrassPoint {
            x,
            y,
            y_infinite: false,
            curve,
        }
    }

    pub fn new_infinite(x: i128, curve: WeierStrass) -> Self {
        WeierStrassPoint {
            x,
            y: i128::MAX,
            y_infinite: true,
            curve,
        }
    }

    pub fn is_infinite(&self) -> bool {
        self.y_infinite
    }

    pub fn print(&self) {
        match self.is_infinite() {
            true => { println!("{}", format!("Point with x={} y=\u{221e}", self.x)) }
            false => { println!("{}", format!("Point with x={} y={}", self.x, self.y)) }
        }
    }

    /// determines the slope of a point and another one
    fn get_slope(&self, other: &WeierStrassPoint) -> Option<i128> {
        // set variables
        let p = self.curve.p;
        let denominator;
        let numerator;

        // determine slope
        if &self == other {
            // point doubling
            denominator = 2 * self.y;
            if denominator == 0 { return None; }

            let modulo = p * denominator;
            numerator = (3 * mod_pow(self.x, 2, modulo) + self.curve.a).rem_euclid(modulo);
        } else {
            // point addition
            denominator = other.x - self.x;
            if denominator == 0 { return None; }

            numerator = (other.y - self.y).rem_euclid(p * denominator)
        }

        // return integer slope
        match mod_inv(denominator, p) {
            Some(inverse) => { Some(mod_mul(numerator, inverse, p)) }
            None => { None }
        }
    }
}


impl PartialEq<WeierStrassPoint> for &WeierStrassPoint {
    fn eq(&self, other: &WeierStrassPoint) -> bool {
        self.x == other.x && self.y == other.y && self.curve == other.curve &&
            self.is_infinite() == other.is_infinite()
    }
}

impl ops::Add<WeierStrassPoint> for WeierStrassPoint {
    type Output = Option<WeierStrassPoint>;

    fn add(self, other: WeierStrassPoint) -> Self::Output {
        // check for matching curves
        if self.curve != other.curve {
            return None;
        }

        // check for infinite points
        if self.is_infinite() {
            return Some(WeierStrassPoint::new_infinite(self.x, self.curve));
        }
        if other.is_infinite() {
            return Some(WeierStrassPoint::new_infinite(self.x, self.curve));
        }

        match self.get_slope(&other) {
            Some(slope) => {
                // determine new coordinates of the new point
                let x = (slope.pow(2) - self.x - other.x).rem_euclid(self.curve.p);
                let y = (slope * (self.x - x) - self.y).rem_euclid(self.curve.p);
                Some(WeierStrassPoint::new(x, y, self.curve))
            }
            None => {
                // infinite slope, so point in infinity is returned
                Some(WeierStrassPoint::new_infinite(other.x, self.curve))
            }
        }
    }
}

impl WeierStrassPoint {
    /// runs one iteration of the lenstra algorithm
    fn lenstra(&self) -> Option<i128> {
        let mut point = self.clone();
        let mut next_point = self.clone();

        let p = self.curve.p;

        // define function
        let mut check_point = |scalar: i128| -> bool {
            let msb_position = get_msb_position(scalar);

            // run a slightly modified version of double and add
            for index in (0..msb_position).rev() {
                // double
                next_point = (point + point).unwrap();
                if next_point.is_infinite() {
                    return true;
                }
                point = next_point;

                // add
                if (scalar >> index) & 0b1 == 0b1 {
                    next_point = (point + self.clone()).unwrap();
                    if next_point.is_infinite() {
                        return true;
                    }
                    point = next_point;
                }
            }
            return false;
        };

        for factorial in 2..=MAX_FACTOR {
            if check_point(factorial) {
                let result = gcd((point.x - next_point.x).rem_euclid(p), p);

                // avoid returning p or 1
                return match p > result && result > 1 {
                    true => { Some(result) }
                    false => { None }
                };
            }
        }

        return None;
    }
}

/// runs the lenstra-factorization algorithm for a provided number
pub fn factorize(number: i128) -> Option<i128> {
    // check for dividable by two
    if (number & 0b1) == 0 {
        return Some(number.div_euclid(2));
    }

    let mut rng = rand::thread_rng();
    for i in 0..MAX_ITERATIONS {
        // get a random curve and point
        let x: i128 = rng.gen_range(0..number.integer_sqrt());
        let y: i128 = rng.gen_range(0..number.integer_sqrt());
        let a: i128 = rng.gen_range(0..number.integer_sqrt());

        let b: i128 = (mod_pow(y, 2, number) - mod_pow(x, 3, number) - a * x).rem_euclid(number);


        let point = match WeierStrass::new(a, b, number) {
            Some(curve) => { WeierStrassPoint::new(x, y, curve) }
            None => { continue; }
        };

        if let Some(factor) = point.lenstra() {
            println!("Finished on {}th iteration", i + 1);
            return Some(factor);
        }
    }

    return None;
}

fn main() {
    let input: i128 = 593 * 1453;

    match factorize(input) {
        Some(result) => { println!("found factor p={}", result) }
        None => { println!("No factors found!") }
    }
}
