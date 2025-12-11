import numpy as np
import matplotlib.pyplot as plt
import sys

# تعریف ثابت اپسیلون برای دقت محاسبات اعداد اعشاری
EPS = sys.float_info.epsilon

# برای مدیریت ذخیره‌سازی و خروجی داده‌های انتگرال‌گیری
class Output:
    
    def __init__(self, nsave=-1):
        # بررسی مقدار nsave برای تعیین حالت ذخیره‌سازی (dense یا غیر dense)
        if nsave == -1:
            self.kmax = -1  # غیرفعال کردن ذخیره‌سازی
            self.dense = False  
        else:
            self.kmax = 500  # حداکثر تعداد نقاط ذخیره‌سازی
            self.nsave = nsave  # تعداد فواصل خروجی
            self.dense = nsave > 0  
        self.count = 0  
        self.xsave = None  
        self.ysave = None 
        self.errsave = None  
        self.hsave = None  

    # مقدار اولیه برای ذخیره
    def init(self, neqn, xlo, xhi):
        self.nvar = neqn  # تعداد متغیرهای معادلات دیفرانسیل
        if self.kmax == -1:
            return  
        self.xsave = np.zeros(self.kmax)  
        self.ysave = np.zeros((self.nvar, self.kmax))  
        self.errsave = np.zeros(self.kmax)  
        self.hsave = np.zeros(self.kmax)  
        if self.dense:
            self.x1 = xlo  # زمان شروع
            self.x2 = xhi  # زمان پایان
            self.xout = self.x1  
            self.dxout = (self.x2 - self.x1) / self.nsave  

    # افزایش اندازه آرایه‌های ذخیره‌سازی در صورت پر شدن
    def resize(self):
        kold = self.kmax  
        self.kmax *= 2  
        tempvec = self.xsave.copy()  
        self.xsave = np.zeros(self.kmax)  
        self.xsave[:kold] = tempvec[:kold]  
        tempmat = self.ysave.copy()  
        self.ysave = np.zeros((self.nvar, self.kmax))  
        self.ysave[:, :kold] = tempmat[:, :kold]  
        temperr = self.errsave.copy()  
        self.errsave = np.zeros(self.kmax)  
        self.errsave[:kold] = temperr[:kold]  
        temph = self.hsave.copy()  
        self.hsave = np.zeros(self.kmax)  
        self.hsave[:kold] = temph[:kold] 

    # ذخیره خروجی متراکم در نقاط مشخص
    def save_dense(self, s, xout, h):
        if self.count == self.kmax:
            self.resize()  # افزایش ظرفیت در صورت پر شدن
        for i in range(self.nvar):
            self.ysave[i, self.count] = s.dense_out(i, xout, h)  # ذخیره مقادیر متغیرها
        self.xsave[self.count] = xout  
        self.errsave[self.count] = s.error()  
        self.hsave[self.count] = h  
        self.count += 1 

    # ذخیره خروجی غیر متراکم
    def save(self, x, y, err, hdid):
        if self.kmax <= 0:
            return  # اگر ذخیره‌سازی غیرفعال باشد، خارج شو
        if self.count == self.kmax:
            self.resize()  # افزایش ظرفیت در صورت پر شدن
        for i in range(self.nvar):
            self.ysave[i, self.count] = y[i]  # ذخیره مقادیر متغیرها
        self.xsave[self.count] = x  
        self.errsave[self.count] = err  
        self.hsave[self.count] = hdid  
        self.count += 1  

    # مدیریت خروجی متراکم در نقاط مورد نظر
    def out(self, nstp, x, y, s, h):
        if not self.dense:
            raise ValueError("dense output not set in Output!")  # خطا در صورت غیرفعال بودن خروجی متراکم
        if nstp == -1:
            self.save(x, y, s.error(), h)  
            self.xout += self.dxout 
        else:
            while (x - self.xout) * (self.x2 - self.x1) > 0.0:
                self.save_dense(s, self.xout, h)  
                self.xout += self.dxout 

# متدهای گام‌زنی عددی
class StepperBase:
    def __init__(self, y, dydx, x, atol, rtol, dense):
        self.n = len(y)  # تعداد معادلات
        self.neqn = self.n  
        self.y = y  
        self.dydx = dydx  
        self.x = x  
        self.atol = atol  # خطای مطلق مجاز
        self.rtol = rtol  # خطای نسبی مجاز
        self.dense = dense  # فعال بودن خروجی متراکم
        self.yout = np.zeros(self.n)  
        self.yerr = np.zeros(self.n)  # خطای محاسبه‌شده
        self.xold = None  
        self.hdid = None  # اندازه گام انجام‌شده
        self.hnext = None  # اندازه گام بعدی

#  رونگ-کوتا دورمند-پرینس مرتبه 5
class StepperDopr5(StepperBase):
    # تنظیم پارامترهای خاص متد دورمند-پرینس
    def __init__(self, y, dydx, x, atol, rtol, dense):
        super().__init__(y, dydx, x, atol, rtol, dense)  
        self.EPS = EPS  
        self.k2 = np.zeros(self.n)  
        self.k3 = np.zeros(self.n)  
        self.k4 = np.zeros(self.n)  
        self.k5 = np.zeros(self.n)  
        self.k6 = np.zeros(self.n) 
        self.rcont1 = np.zeros(self.n)  
        self.rcont2 = np.zeros(self.n)  
        self.rcont3 = np.zeros(self.n)  
        self.rcont4 = np.zeros(self.n)  
        self.rcont5 = np.zeros(self.n)  
        self.dydxnew = np.zeros(self.n)  
        self.con = self.Controller()  

    # کنترل اندازه گام
    class Controller:
        # متد سازنده برای کنترلر
        def __init__(self):
            self.reject = False  # وضعیت رد گام
            self.errold = 1.0e-4  # خطای قبلی
            self.hnext = None  # اندازه گام بعدی

        # بررسی موفقیت گام و تنظیم اندازه گام بعدی
        def success(self, err, h):
            beta = 0.0  # پارامتر بتا برای کنترل خطا
            alpha = 0.2 - beta * 0.75  # پارامتر آلفا
            safe = 0.9  # ضریب ایمنی
            minscale = 0.2  # حداقل مقیاس گام
            maxscale = 10.0  # حداکثر مقیاس گام
            if err <= 1.0:
                if err == 0.0:
                    scale = maxscale  # حداکثر مقیاس در صورت خطای صفر
                else:
                    scale = safe * err ** (-alpha) * self.errold ** beta  # محاسبه مقیاس
                    if scale < minscale:
                        scale = minscale  # محدود کردن به حداقل
                    if scale > maxscale:
                        scale = maxscale  # محدود کردن به حداکثر
                if self.reject:
                    self.hnext = h * min(scale, 1.0)  # محدود کردن گام در صورت رد قبلی
                else:
                    self.hnext = h * scale  # تنظیم گام بعدی
                self.errold = max(err, 1.0e-4)  # به‌روزرسانی خطای قبلی
                self.reject = False 
                return True, h
            else:
                scale = max(safe * err ** (-alpha), minscale)  # کاهش گام در صورت خطای زیاد
                new_h = h * scale  # محاسبه گام جدید
                self.reject = True  # علامت‌گذاری گام به‌عنوان رد شده
                return False, new_h

    #  رونگ-کوتا مرتبه 5
    def dy(self, h, derivs):
        # ضرایب روش دورمند-پرینس
        c2 = 0.2
        c3 = 0.3
        c4 = 0.8
        c5 = 8.0 / 9.0
        a21 = 0.2
        a31 = 3.0 / 40.0
        a32 = 9.0 / 40.0
        a41 = 44.0 / 45.0
        a42 = -56.0 / 15.0
        a43 = 32.0 / 9.0
        a51 = 19372.0 / 6561.0
        a52 = -25360.0 / 2187.0
        a53 = 64448.0 / 6561.0
        a54 = -212.0 / 729.0
        a61 = 9017.0 / 3168.0
        a62 = -355.0 / 33.0
        a63 = 46732.0 / 5247.0
        a64 = 49.0 / 176.0
        a65 = -5103.0 / 18656.0
        a71 = 35.0 / 384.0
        a73 = 500.0 / 1113.0
        a74 = 125.0 / 192.0
        a75 = -2187.0 / 6784.0
        a76 = 11.0 / 84.0
        e1 = 71.0 / 57600.0
        e3 = -71.0 / 16695.0
        e4 = 71.0 / 1920.0
        e5 = -17253.0 / 339200.0
        e6 = 22.0 / 525.0
        e7 = -1.0 / 40.0

       
        ytemp = self.y + h * a21 * self.dydx
        self.k2 = derivs(self.x + c2 * h, ytemp)

        ytemp = self.y + h * (a31 * self.dydx + a32 * self.k2)
        self.k3 = derivs(self.x + c3 * h, ytemp)

        ytemp = self.y + h * (a41 * self.dydx + a42 * self.k2 + a43 * self.k3)
        self.k4 = derivs(self.x + c4 * h, ytemp)

        ytemp = self.y + h * (a51 * self.dydx + a52 * self.k2 + a53 * self.k3 + a54 * self.k4)
        self.k5 = derivs(self.x + c5 * h, ytemp)

        ytemp = self.y + h * (a61 * self.dydx + a62 * self.k2 + a63 * self.k3 + a64 * self.k4 + a65 * self.k5)
        xph = self.x + h
        self.k6 = derivs(xph, ytemp)

        # محاسبه خروجی و خطا
        self.yout = self.y + h * (a71 * self.dydx + a73 * self.k3 + a74 * self.k4 + a75 * self.k5 + a76 * self.k6)
        self.dydxnew = derivs(xph, self.yout)
        self.yerr = h * (e1 * self.dydx + e3 * self.k3 + e4 * self.k4 + e5 * self.k5 + e6 * self.k6 + e7 * self.dydxnew)

    # محاسبه خطای گام
    def error(self):
        err = 0.0
        for i in range(self.n):
            sk = self.atol + self.rtol * max(abs(self.y[i]), abs(self.yout[i]))  # مقیاس خطا
            err += (self.yerr[i] / sk) ** 2  # جمع مربعات خطاها
        return np.sqrt(err / self.n)  # خطای نرمال‌شده

    # آماده‌سازی خروجی متراکم
    def prepare_dense(self, h, derivs):
        # تعریف ضرایب برای درون‌یابی متراکم
        d1 = -12715105075.0 / 11282082432.0
        d3 = 87487479700.0 / 32700410799.0
        d4 = -10690763975.0 / 1880347072.0
        d5 = 701980252875.0 / 199316789632.0
        d6 = -1453857185.0 / 822651844.0
        d7 = 69997945.0 / 29380423.0
        self.rcont1 = self.y.copy()  
        ydiff = self.yout - self.y  
        self.rcont2 = ydiff.copy()
        bspl = h * self.dydx - ydiff
        self.rcont3 = bspl.copy()
        self.rcont4 = ydiff - h * self.dydxnew - bspl
        self.rcont5 = h * (d1 * self.dydx + d3 * self.k3 + d4 * self.k4 + d5 * self.k5 + d6 * self.k6 + d7 * self.dydxnew)

    # محاسبه مقدار خروجی متراکم در نقطه خاص
    def dense_out(self, i, x, h):
        s = (x - self.xold) / h  # پارامتر درون‌یابی
        s1 = 1.0 - s
        return self.rcont1[i] + s * (self.rcont2[i] + s1 * (self.rcont3[i] + s * (self.rcont4[i] + s1 * self.rcont5[i])))

    # قدم بردار
    def step(self, htry, derivs):
        h = htry  # اندازه گام پیشنهادی
        for _ in range(100):
            self.dy(h, derivs)  # محاسبه گام
            err = self.error()  # محاسبه خطا
            success, h = self.con.success(err, h)  # بررسی موفقیت گام
            if success:
                break
            if abs(h) <= abs(self.x) * self.EPS:
                raise ValueError("stepsize underflow in StepperDopr5")  # خطا در صورت کوچک شدن بیش از حد گام
        if self.dense:
            self.prepare_dense(h, derivs)  
        self.dydx = self.dydxnew.copy()  
        self.y = self.yout.copy()  
        self.xold = self.x  
        self.x += h  
        self.hdid = h  
        self.hnext = self.con.hnext 

# کلاس اصلی برای انتگرال‌گیری معادلات دیفرانسیل
class Odeint:
    MAXSTP = 50000  # حداکثر تعداد گام‌ها

    # متد سازنده برای تنظیم پارامترهای انتگرال‌گیری
    def __init__(self, ystart, x1, x2, atol, rtol, h1, hmin, out, derivs):
        self.EPS = EPS  # دقت اعشاری سیستم
        self.nok = 0  # گام خوب
        self.nbad = 0  # گام بد
        self.nvar = len(ystart)  
        self.x1 = x1  
        self.x2 = x2  
        self.hmin = hmin  # حداقل اندازه گام
        self.dense = out.dense  
        self.ystart = ystart  # شرایط اولیه
        self.out = out  
        self.derivs = derivs  # تابع معادلات دیفرانسیل
        self.y = np.array(ystart)  
        self.dydx = np.zeros(self.nvar)  
        self.x = x1  
        direction = self.x2 - self.x1  
        self.h = h1 if direction == 0 else np.sign(h1) * abs(h1) if direction * h1 >= 0 else -abs(h1)  # تنظیم اندازه گام اولیه
        self.s = StepperDopr5(self.y, self.dydx, self.x, atol, rtol, self.dense)  
        self.out.init(self.s.neqn, self.x1, self.x2) 

    # انتگرال‌گیری
    def integrate(self):
        self.dydx = self.derivs(self.x, self.y)  # محاسبه مشتقات اولیه
        if self.dense:
            self.out.out(-1, self.x, self.y, self.s, self.h)  # ذخیره خروجی اولیه در حالت متراکم
        else:
            self.out.save(self.x, self.y, self.s.error(), self.s.hdid)  # ذخیره خروجی اولیه
        self.nstp = 0  # شمارشگر گام‌ها
        while self.nstp < self.MAXSTP:
            if (self.x + self.h * 1.0001 - self.x2) * (self.x2 - self.x1) > 0.0:
                self.h = self.x2 - self.x  # تنظیم گام برای رسیدن به نقطه پایان
            self.s.step(self.h, self.derivs)  # اجرای گام
            self.x = self.s.x  
            self.y = self.s.y  
            if self.s.hdid == self.h:
                self.nok += 1  # گام خوب
            else:
                self.nbad += 1  # گام بد
            if self.dense:
                self.out.out(self.nstp, self.x, self.y, self.s, self.s.hdid)  
            else:
                self.out.save(self.x, self.y, self.s.error(), self.s.hdid)  
            if (self.x - self.x2) * (self.x2 - self.x1) >= 0.0:
                self.ystart[:] = self.y  
                if self.out.kmax > 0 and abs(self.out.xsave[self.out.count - 1] - self.x2) > 100.0 * abs(self.x2) * self.EPS:
                    self.out.save(self.x, self.y, self.s.error(), self.s.hdid)  
                return
            if abs(self.s.hnext) <= self.hmin:
                raise ValueError("Step size too small in Odeint")  # خطا در صورت کوچک شدن بیش از حد گام
            self.h = self.s.hnext  # به‌روزرسانی اندازه گام
            self.nstp += 1  # افزایش شمارشگر گام‌ها
        raise ValueError("Too many steps in routine Odeint")  # خطا در صورت تعداد گام‌های زیاد

# تعریف تابع معادلات دیفرانسیل برای نوسانگر دافینگ
def f(t, y):
    x, v, s = y  
    dxdt = v  
    dvdt = -0.2 * v + x - 0.1 * (x ** 3) + np.cos(t)  
    dsdt = np.cos(t / 2)  
    return np.array([dxdt, dvdt, dsdt])  

# ثابت ها
x1 = 0.0  
x2 = 400.0  
ystart = np.array([0.0, 0.0, 0.0])  # شرایط اولیه
atol = 1e-4  # خطای مطلق
rtol = 1e-4  # خطای نسبی
h1 = 0.1  # اندازه گام اولیه
hmin = 0.0  # حداقل اندازه گام
nsave = 4000  # تعداد فواصل خروجی متراکم

# انتگرال‌گیری
out = Output(nsave)
ode = Odeint(ystart, x1, x2, atol, rtol, h1, hmin, out, f)
ode.integrate()  # اجرای انتگرال‌گیری

# استخراج نتایج
t = out.xsave[:out.count]  
x = out.ysave[0, :out.count]  
v = out.ysave[1, :out.count]  
s = out.ysave[2, :out.count]  
err = out.errsave[:out.count]  
h = out.hsave[:out.count]  

# محاسبه بخش پوآنکاره در نقاط t = 2pi n
poincare_x = []
poincare_v = []
t_poincare = np.arange(0, x2, 2 * np.pi)  # زمان‌های t = 2pi n
for tp in t_poincare[1:]:  # نادیده گرفتن t=0 برای اجتناب از اثرات گذرا
    idx = np.argmin(np.abs(t - tp))  # یافتن نزدیک‌ترین شاخص
    if abs(t[idx] - tp) < 0.01:  # اطمینان از نزدیکی به زمان مورد نظر
        poincare_x.append(x[idx]) 
        poincare_v.append(v[idx])  

# رسم نتایج 
fig, axs = plt.subplots(3, 2, figsize=(12, 12))  

axs[0, 0].plot(t, x)
axs[0, 0].set_title('x(t)')
axs[0, 0].set_xlabel('t')
axs[0, 0].set_ylabel('x')
axs[0, 0].grid(True, linestyle='--', alpha=0.7)

axs[0, 1].plot(x, v)
axs[0, 1].set_title('Phase Space (x vs v)')
axs[0, 1].set_xlabel('x')
axs[0, 1].set_ylabel('v')
axs[0, 1].grid(True, linestyle='--', alpha=0.7)

axs[1, 0].semilogy(t, err, marker='', linestyle='-')
axs[1, 0].set_title('Error in each step')
axs[1, 0].set_xlabel('t')
axs[1, 0].set_ylabel('errmax')
axs[1, 0].grid(True, linestyle='--', alpha=0.7)

axs[1, 1].semilogy(t, h, marker='', linestyle='-')
axs[1, 1].set_title('Step Size')
axs[1, 1].set_xlabel('t')
axs[1, 1].set_ylabel('Step Size')
axs[1, 1].grid(True, linestyle='--', alpha=0.7)

axs[2, 0].scatter(poincare_x, poincare_v, s=10, c='red')
axs[2, 0].set_title('Poincaré Section')
axs[2, 0].set_xlabel('x')
axs[2, 0].set_ylabel('v')
axs[2, 0].grid(True, linestyle='--', alpha=0.7)

axs[2, 1] = fig.add_subplot(3, 2, 6, projection='3d')  
axs[2, 1].plot(t, x, v)  
axs[2, 1].set_title('3D Trajectory (t, x, v)')  
axs[2, 1].set_xlabel('t')  
axs[2, 1].set_ylabel('x')  
axs[2, 1].set_zlabel('v')  
axs[2, 1].grid(True, linestyle='--', alpha=0.7) 

plt.tight_layout()
plt.show()