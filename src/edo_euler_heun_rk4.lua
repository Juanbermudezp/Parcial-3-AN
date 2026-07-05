platform.apiLevel = '2.0'

-- Programa TI-Nspire: solucion numerica de EDO de primer orden
-- y' = f(x,y), y(x0) = y0, por los metodos de Euler, Heun y Runge-Kutta 4

local MAX_N = 150
local ROW_H = 14
local SCROLL_STEP = 3

local fields = {
    { label = "f(x,y)=", value = "x+y" },
    { label = "x0=",     value = "0" },
    { label = "y0=",     value = "1" },
    { label = "xn=",     value = "1" },
    { label = "h=",      value = "0.1" },
}
local focus = 1

local results, errorMsg, scrollY = nil, nil, 0

local function toNumber(s)
    return tonumber(s)
end

local function fmt(v)
    return string.format("%.6g", v)
end

-- Evalua f(x,y) en el servidor de matematicas (CAS) sustituyendo x,y por
-- valores numericos en notacion fija (evita el sufijo cientifico "e" que el
-- CAS interpretaria como el numero de Euler en lugar de exponente).
local function feval(expr, xv, yv)
    local s = string.format("(%s)|x=%.12f and y=%.12f", expr, xv, yv)
    local res, err = math.eval(s)
    if res == nil then
        return nil, "No se pudo evaluar f(x,y)"
    end
    if type(res) ~= "number" then
        return nil, "f(x,y) debe producir un numero"
    end
    return res
end

local function euler(expr, x0, y0, h, n)
    local xs, ys = { [0] = x0 }, { [0] = y0 }
    for i = 0, n - 1 do
        local f1, err = feval(expr, xs[i], ys[i])
        if not f1 then return nil, nil, err end
        xs[i + 1] = xs[i] + h
        ys[i + 1] = ys[i] + h * f1
    end
    return xs, ys
end

local function heun(expr, x0, y0, h, n)
    local xs, ys = { [0] = x0 }, { [0] = y0 }
    for i = 0, n - 1 do
        local f1, e1 = feval(expr, xs[i], ys[i])
        if not f1 then return nil, nil, e1 end
        local xnext = xs[i] + h
        local ypred = ys[i] + h * f1
        local f2, e2 = feval(expr, xnext, ypred)
        if not f2 then return nil, nil, e2 end
        xs[i + 1] = xnext
        ys[i + 1] = ys[i] + (h / 2) * (f1 + f2)
    end
    return xs, ys
end

local function rk4(expr, x0, y0, h, n)
    local xs, ys = { [0] = x0 }, { [0] = y0 }
    for i = 0, n - 1 do
        local k1, e1 = feval(expr, xs[i], ys[i])
        if not k1 then return nil, nil, e1 end
        local k2, e2 = feval(expr, xs[i] + h / 2, ys[i] + h / 2 * k1)
        if not k2 then return nil, nil, e2 end
        local k3, e3 = feval(expr, xs[i] + h / 2, ys[i] + h / 2 * k2)
        if not k3 then return nil, nil, e3 end
        local k4, e4 = feval(expr, xs[i] + h, ys[i] + h * k3)
        if not k4 then return nil, nil, e4 end
        xs[i + 1] = xs[i] + h
        ys[i + 1] = ys[i] + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    end
    return xs, ys
end

local function compute()
    errorMsg, results, scrollY = nil, nil, 0

    local expr = fields[1].value
    local x0 = toNumber(fields[2].value)
    local y0 = toNumber(fields[3].value)
    local xn = toNumber(fields[4].value)
    local h = toNumber(fields[5].value)

    if expr == "" or not x0 or not y0 or not xn or not h then
        errorMsg = "Verifique los datos ingresados"
        return
    end
    if h == 0 then
        errorMsg = "h no puede ser 0"
        return
    end
    if (xn - x0) / h < 0 then
        errorMsg = "El signo de h no coincide con [x0,xn]"
        return
    end

    local n = math.floor((xn - x0) / h + 0.5)
    if n < 1 then
        errorMsg = "Rango invalido"
        return
    end
    if n > MAX_N then
        errorMsg = "Demasiados pasos (maximo " .. MAX_N .. ")"
        return
    end

    local xe, ye, errE = euler(expr, x0, y0, h, n)
    if not xe then errorMsg = errE; return end
    local xh, yh, errH = heun(expr, x0, y0, h, n)
    if not xh then errorMsg = errH; return end
    local xr, yr, errR = rk4(expr, x0, y0, h, n)
    if not xr then errorMsg = errR; return end

    results = {}
    for i = 0, n do
        results[#results + 1] = { x = xe[i], e = ye[i], h = yh[i], r = yr[i] }
    end
end

function on.paint(gc)
    local w, hgt = platform.window:width(), platform.window:height()
    gc:setColorRGB(255, 255, 255)
    gc:fillRect(0, 0, w, hgt)
    gc:setColorRGB(0, 0, 0)

    gc:setFont("sansserif", "b", 10)
    gc:drawString("EDO: Euler / Heun / Runge-Kutta 4", 2, 2, "top")

    gc:setFont("sansserif", "r", 9)
    local y = 18
    local labelW = 42
    for i, f in ipairs(fields) do
        gc:drawString(f.label, 2, y, "top")
        local boxX, boxW = labelW, w - labelW - 4
        gc:setColorRGB(focus == i and 205 or 235, focus == i and 225 or 235, focus == i and 255 or 235)
        gc:fillRect(boxX, y, boxW, ROW_H)
        gc:setColorRGB(0, 0, 0)
        gc:drawRect(boxX, y, boxW, ROW_H)
        gc:drawString(f.value, boxX + 2, y, "top")
        y = y + ROW_H + 2
    end

    gc:setFont("sansserif", "r", 7)
    gc:drawString("Enter=calcular  Tab=cambiar campo  Arriba/Abajo=desplazar tabla", 2, y, "top")
    y = y + 12

    gc:setFont("sansserif", "r", 9)
    if errorMsg then
        gc:setColorRGB(200, 0, 0)
        gc:drawString(errorMsg, 2, y, "top")
        gc:setColorRGB(0, 0, 0)
        return
    end

    if not results then
        gc:drawString("Ingrese los datos y presione Enter", 2, y, "top")
        return
    end

    local cols = {
        { name = "i",     w = 20 },
        { name = "x",     w = 55 },
        { name = "Euler", w = 70 },
        { name = "Heun",  w = 70 },
        { name = "RK4",   w = 70 },
    }

    gc:setColorRGB(220, 220, 220)
    gc:fillRect(2, y, w - 4, ROW_H)
    gc:setColorRGB(0, 0, 0)
    local cx = 2
    for _, c in ipairs(cols) do
        gc:drawString(c.name, cx, y, "top")
        cx = cx + c.w
    end
    y = y + ROW_H
    gc:drawLine(2, y, w - 2, y)

    local visibleRows = math.floor((hgt - y) / ROW_H)
    local lastRow = math.min(#results, scrollY + visibleRows)
    for r = scrollY + 1, lastRow do
        local row = results[r]
        local vals = { tostring(r - 1), fmt(row.x), fmt(row.e), fmt(row.h), fmt(row.r) }
        cx = 2
        for ci, c in ipairs(cols) do
            gc:drawString(vals[ci], cx, y, "top")
            cx = cx + c.w
        end
        y = y + ROW_H
    end
end

function on.charIn(ch)
    local f = fields[focus]
    f.value = f.value .. ch
    platform.window:invalidate()
end

function on.backspaceKey()
    local f = fields[focus]
    f.value = f.value:sub(1, -2)
    platform.window:invalidate()
end

function on.tabKey()
    focus = focus % #fields + 1
    platform.window:invalidate()
end

function on.enterKey()
    compute()
    platform.window:invalidate()
end

function on.arrowKey(key)
    if not results then return end
    local maxScroll = math.max(0, #results - 1)
    if key == "down" then
        scrollY = math.min(scrollY + SCROLL_STEP, maxScroll)
    elseif key == "up" then
        scrollY = math.max(scrollY - SCROLL_STEP, 0)
    end
    platform.window:invalidate()
end

function on.mouseDown(x, y)
    local yy = 18
    for i in ipairs(fields) do
        if y >= yy and y <= yy + ROW_H then
            focus = i
            platform.window:invalidate()
            return
        end
        yy = yy + ROW_H + 2
    end
end

function on.save()
    local t = {}
    for i, f in ipairs(fields) do t[i] = f.value end
    return t
end

function on.restore(t)
    if not t then return end
    for i, f in ipairs(fields) do
        if t[i] then f.value = t[i] end
    end
end
