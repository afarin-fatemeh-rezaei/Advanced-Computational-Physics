import numpy as np
import matplotlib.pyplot as plt
import sys


# هدایت خروجی به فایل
sys.stdout = open('main.txt', 'w')
sys.stderr = sys.stdout

# ثابت‌های مورد استفاده در الگوریتم
SAFETY = 0.9  # ضریب ایمنی
PGRDW = -0.2 # توان بزرگ کننده
PSHRNK = -0.25 # توان کوچک کننده
ERRCON = (5.0 / SAFETY) ** (1.0 / PGRDW)
MAXSTP = 10000 # حداکثر تعداد گام‌ها
TINY = 1.0e-30 # عدد ایمنی 

# تابع معادلات دیفرانسیل برای نوسانگر دافینگ با متغیر کمکی
def f(t, y):
    x, v, s = y  
    dxdt = v
    dvdt = -(0.2 * v) + x - (0.1 * (x**3)) + np.cos(t)
    dsdt = np.cos(t / 2)  # مشتق متغیر کمکی
    return np.array([dxdt, dvdt, dsdt])

# روش رونگ-کوتا مرتبه پنجم
def rkck(y, dydx, n, x, h, derivs):
    # ضرایب کش-کارپ
    a2, a3, a4, a5, a6 = 0.2, 0.3, 0.6, 1.0, 0.875
    b21 = 0.2
    b31, b32 = 3.0/40.0, 9.0/40.0
    b41, b42, b43 = 0.3, -0.9, 1.2
    b51, b52, b53, b54 = -11.0/54.0, 2.5, -70.0/27.0, 35.0/27.0
    b61, b62, b63, b64, b65 = 1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0
    c1, c3, c4, c6 = 37.0/378.0, 250.0/621.0, 125.0/594.0, 512.0/1771.0
    dc1, dc3, dc4, dc5, dc6 = c1 - 2825.0/27648.0, c3 - 18575.0/48384.0, c4 - 13525.0/55296.0, -277.0/14336.0, c6 - 0.25
    
    ak2 = np.zeros(n)
    ak3 = np.zeros(n)
    ak4 = np.zeros(n)
    ak5 = np.zeros(n)
    ak6 = np.zeros(n)
    ytemp = np.zeros(n)
    
    for i in range(n):
        ytemp[i] = y[i] + b21 * h * dydx[i]
    ak2 = derivs(x + a2 * h, ytemp)
    
    for i in range(n):
        ytemp[i] = y[i] + h * (b31 * dydx[i] + b32 * ak2[i])
    ak3 = derivs(x + a3 * h, ytemp)
    
    for i in range(n):
        ytemp[i] = y[i] + h * (b41 * dydx[i] + b42 * ak2[i] + b43 * ak3[i])
    ak4 = derivs(x + a4 * h, ytemp)
    
    for i in range(n):
        ytemp[i] = y[i] + h * (b51 * dydx[i] + b52 * ak2[i] + b53 * ak3[i] + b54 * ak4[i])
    ak5 = derivs(x + a5 * h, ytemp)
    
    for i in range(n):
        ytemp[i] = y[i] + h * (b61 * dydx[i] + b62 * ak2[i] + b63 * ak3[i] + b64 * ak4[i] + b65 * ak5[i])
    ak6 = derivs(x + a6 * h, ytemp)
    
    yout = np.zeros(n)
    for i in range(n):
        yout[i] = y[i] + h * (c1 * dydx[i] + c3 * ak3[i] + c4 * ak4[i] + c6 * ak6[i])
    
    yerr = np.zeros(n)
    for i in range(n):
        yerr[i] = h * (dc1 * dydx[i] + dc3 * ak3[i] + dc4 * ak4[i] + dc5 * ak5[i] + dc6 * ak6[i])

    print(f"\n\nrkck1. yerr_x={abs(yerr[0])}, yerr_v={abs(yerr[1])}, yerr_s={abs(yerr[2])}, yout={yout}")
    return yout, yerr



# تابع تنظیم گام برای کنترل خطا
def rkqs(y, dydx, n, x, htry, eps, yscal, derivs):
    h = htry
    while True:
        ytemp, yerr = rkck(y, dydx, n, x, h, derivs)
        print(f"\nrkqs1. h={h:.6e}")

        errmax = 0
        for i in range(n):
            errmax = max(errmax, np.max(np.abs(yerr) /  (np.abs(yscal))))

        print(f"\nrkqs2. errmax={abs(errmax)}, yerr/yscal={np.max(np.abs(yerr)/ np.abs(yscal))}, yscal={yscal}")
        
        errmax /= eps
        print(f"\nrkqs3. errmax={abs(errmax)}")
        if errmax <= 1:
            break
        
        htemp = SAFETY * h * errmax ** PSHRNK # کوچک کننده گام
        h = htemp
        print(f"\nrkqs4. htemp={htemp}")
        if h >= 0.0:
            h = max(htemp, 0.1 * h) # گام تا یک حد کوچک تر نشود 
        else:
            h = min(htemp, 0.1 * h)
        print(f"\nrkqs5. h={h}")
        
        xnew = x + h
        if xnew == x:
            raise RuntimeError("stepsize underflow in rkqs")
    
    if errmax > ERRCON:
        hnext = SAFETY * h * errmax ** PGRDW #بزرگ کننده گام
    else:
        hnext = 5 * h  # گام تا یک حد بزرگ تر نشود

    print(f"\nrkqs6. hnext={hnext}")
    return ytemp, x + h, h, hnext, errmax

# تابع برای تعیین علامت مقدار
def SIGN(a, b):
    return abs(a) if b >= 0 else -abs(a)

# تابع اصلی
def odeint(ystart, x1, x2, eps, h1, hmin, derivs, rkqs_func, dxsav, kmax):
    nvar = len(ystart)
    y = np.array(ystart, dtype=float)   # کپی از شرایط اولیه
    x = float(x1)
    h = SIGN(h1, x2 - x1)
    nok = 0
    nbad = 0
    kount = 0
    xp = []  
    yp = [[] for _ in range(nvar)]  
    xsav = x - dxsav * 2.0 if kmax > 0 else x
    err_max = []
    steps = []
    errmax = 0

    # ذخیره‌سازی اولیه نتایج
    if kmax > 0:
        xp.append(x)
        for i in range(nvar):
            yp[i].append(y[i])
        err_max.append(0.0)
        steps.append(0.0)

    for nstp in range(MAXSTP):
        dydx = derivs(x, y)
        yscal = np.zeros(nvar)
        for i in range(nvar):
            yscal[i] = abs(y[i]) + abs(dydx[i] * h) + TINY

        # ذخیره نتایج میانی
        if kmax > 0 and kount < kmax - 1 and abs(x - xsav) > abs(dxsav):
            for i in range(nvar):
                yp[i].append(y[i])
            xp.append(x)
            err_max.append(errmax)
            steps.append(h)
            kount += 1
            xsav = x

        # تنظیم اندازه گام در صورت بزرگ بودن زیاد
        if (x + h - x2) * (x + h - x1) > 0.0:
            h = x2 - x

        # اجرای رونگ-کوتا
        y, x_new, hdid, hnext, errmax = rkqs_func(y, dydx, nvar, x, h, eps, yscal, derivs)
        x = x_new

        # تعداد گام های خوب و بد
        if hdid == h:
            nok += 1
        else:
            nbad += 1

        print(f"\nodeint: step {nstp}, x={x:.6f}, hnext={hnext:.2e}, nok={nok}, nbad={nbad}, y={y}")

        # بررسی اتمام انتگرال‌گیری
        if (x - x2) * (x2 - x1) >= 0.0:
            for i in range(nvar):
                ystart[i] = y[i]  
            if kmax != 0:
                for i in range(nvar):
                    yp[i].append(y[i])
                xp.append(x)
                err_max.append(errmax)
                steps.append(h)
                kount += 1
            return np.array(xp), np.array(yp).T, nok, nbad, np.array(err_max), np.array(steps)

        # بررسی اندازه گام نسبت به حداقل
        if abs(hnext) <= hmin:
            raise RuntimeError("Step size too small in odeint")
            
        h = hnext

    raise RuntimeError("Too many steps in routine odeint")

# پارامترهای اولیه
t0 = 0.0
t_end = 400.0 # غیر قابل تغییر 
dt_store = 0.1 # غیر قابل تغییر
eps = 1e-4 # غیر قابل تغییر 
h1 = 0.1 # اندازه گام اولیه
y0 = np.array([0.0, 0.0, 0.0])  # شرایط اولیه غیر قابل تغییر
hmin = 0 # حداقل اندازه گام
kmax = 100000 # بیشترین مقدار ذخیره

# اجرای انتگرال‌گیری
t_vals, y_vals, nok, nbad, err_max, steps = odeint(y0, t0, t_end, eps, h1, hmin, f, rkqs, dt_store, kmax)

# # محاسبه بخش پوآنکاره در نقاط t = 2pi n
poincare_x = []
poincare_v = []
tolerance = 0.1  # تلورانس تشخیص 

for i in range(1, len(t_vals)):
    s_current = y_vals[i, 2]  
    s_previous = y_vals[i-1, 2]  
    t_current = t_vals[i] 
    # # بررسی عبور از صفر و شرط ds/dt > 0 
    if s_previous * s_current <= 0 and np.cos(t_current / 2) > 0: # زمان‌های t = 2pi n
        # درون‌یابی برای یافتن نقطه صفر دقیق
        t_interp = t_vals[i-1] + (t_vals[i] - t_vals[i-1]) * (-s_previous) / (s_current - s_previous)
        x_interp = y_vals[i-1, 0] + (y_vals[i, 0] - y_vals[i-1, 0]) * (-s_previous) / (s_current - s_previous)
        v_interp = y_vals[i-1, 1] + (y_vals[i, 1] - y_vals[i-1, 1]) * (-s_previous) / (s_current - s_previous)
        poincare_x.append(x_interp)
        poincare_v.append(v_interp)

# رسم نتایج
fig, axs = plt.subplots(3, 2, figsize=(12, 12))  

axs[0, 0].plot(t_vals, y_vals[:, 0])
axs[0, 0].set_title('x(t)')
axs[0, 0].set_xlabel('t')
axs[0, 0].set_ylabel('x')
axs[0, 0].grid(True, linestyle='--', alpha=0.7)

axs[0, 1].plot(y_vals[:, 0], y_vals[:, 1])
axs[0, 1].set_title('Phase Space (x vs v)')
axs[0, 1].set_xlabel('x')
axs[0, 1].set_ylabel('v')
axs[0, 1].grid(True, linestyle='--', alpha=0.7)

axs[1, 0].semilogy(t_vals, err_max, marker='', linestyle='-')
axs[1, 0].set_title('Error (errmax)')
axs[1, 0].set_xlabel('t')
axs[1, 0].set_ylabel('errmax')
axs[1, 0].grid(True, linestyle='--', alpha=0.7)

axs[1, 1].semilogy(t_vals, steps, marker='', linestyle='-')
axs[1, 1].set_title('Step Size')
axs[1, 1].set_xlabel('t')
axs[1, 1].set_ylabel('Step Size')
axs[1, 1].grid(True, linestyle='--', alpha=0.7)

axs[2, 0].scatter(poincare_x, poincare_v, s=10, c='red')
axs[2, 0].set_title('Poincaré Section')
axs[2, 0].set_xlabel('x')
axs[2, 0].set_ylabel('v')
axs[2, 0].grid(True, linestyle='--', alpha=0.7)
axs[2, 0].legend()

axs[2, 1] = fig.add_subplot(3, 2, 6, projection='3d')  
axs[2, 1].plot(t_vals, y_vals[:, 0], y_vals[:, 1])  
axs[2, 1].set_title('3D Trajectory (t, x, v)')  
axs[2, 1].set_xlabel('t')  
axs[2, 1].set_ylabel('x')  
axs[2, 1].set_zlabel('v')  
axs[2, 1].grid(True, linestyle='--', alpha=0.7) 

plt.tight_layout()
plt.show()

# چاپ نتایج نهایی
print(f"\nIntegration completed!")
print(f"Steps taken: {len(t_vals)}")
print(f"Successful steps: {nok}, Failed steps: {nbad}")
print(f"Final position: {y_vals[-1, 0]:.6f}, Final velocity: {y_vals[-1, 1]:.6f}, Final s: {y_vals[-1, 2]:.6f}")
print(f"Poincaré points: {len(poincare_x)}")
print("done!")