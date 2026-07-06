platform.apiLevel = '2.0'

-- Programa TI-Nspire: solucion numerica de EDO de primer orden
-- y' = f(x,y), y(x0) = y0, por el metodo multipaso de
-- Adams-Bashforth-Moulton de 4 pasos (predictor-corrector).
-- Los primeros 4 puntos se generan con Runge-Kutta 4 (arranque).
-- Muestra el desarrollo paso a paso y al final la tabla de resultados.

local MAX_N = 150
local ROW_H = 12
local SCROLL_STEP = 3

local fields = {
    { label = "f(x,y)=", value = "x+y" },
    { label = "x0=",     value = "0" },
    { label = "y0=",     value = "1" },
    { label = "xn=",     value = "1" },
    { label = "h=",      value = "0.1" },
}
local focus = 1

local lines, errorMsg, scrollY = nil, nil, 0

local function toNumber(s)
    return tonumber(s)
end

local function fmt(v)
    return string.format("%.6f", v)
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

local function solve(expr, x0, y0, h, n)
    local xs, ys, fs = { [0] = x0 }, { [0] = y0 }, {}
    local log = { "Arranque con Runge-Kutta 4 (i=0..3):" }

    local f0, e0 = feval(expr, xs[0], ys[0])
    if not f0 then return nil, nil, nil, e0 end
    fs[0] = f0

    for i = 0, 2 do
        local k1, e1 = feval(expr, xs[i], ys[i])
        if not k1 then return nil, nil, nil, e1 end
        local xm, ym1 = xs[i] + h / 2, ys[i] + h / 2 * k1
        local k2, e2 = feval(expr, xm, ym1)
        if not k2 then return nil, nil, nil, e2 end
        local ym2 = ys[i] + h / 2 * k2
        local k3, e3 = feval(expr, xm, ym2)
        if not k3 then return nil, nil, nil, e3 end
        local xnext, yp4 = xs[i] + h, ys[i] + h * k3
        local k4, e4 = feval(expr, xnext, yp4)
        if not k4 then return nil, nil, nil, e4 end
        local ynext = ys[i] + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

        log[#log + 1] = string.format("  RK4 paso %d: k1=%.6f k2=%.6f k3=%.6f k4=%.6f", i + 1, k1, k2, k3, k4)
        log[#log + 1] = string.format("  y%d = %.6f   (x%d=%.4f)", i + 1, ynext, i + 1, xnext)

        xs[i + 1], ys[i + 1] = xnext, ynext
        local fnext, efn = feval(expr, xnext, ynext)
        if not fnext then return nil, nil, nil, efn end
        fs[i + 1] = fnext
    end

    log[#log + 1] = ""
    log[#log + 1] = "Adams-Bashforth-Moulton (predictor-corrector):"

    for i = 3, n - 1 do
        local yp = ys[i] + (h / 24) * (55 * fs[i] - 59 * fs[i - 1] + 37 * fs[i - 2] - 9 * fs[i - 3])
        local xnext = xs[i] + h
        local fp, efp = feval(expr, xnext, yp)
        if not fp then return nil, nil, nil, efp end
        local ynext = ys[i] + (h / 24) * (9 * fp + 19 * fs[i] - 5 * fs[i - 1] + fs[i - 2])
        local fnext, efn = feval(expr, xnext, ynext)
        if not fnext then return nil, nil, nil, efn end

        log[#log + 1] = string.format("Paso %d:", i + 1)
        log[#log + 1] = string.format("  predictor yp=y%d+(h/24)(55f%d-59f%d+37f%d-9f%d)", i, i, i - 1, i - 2, i - 3)
        log[#log + 1] = string.format("  yp = %.6f", yp)
        log[#log + 1] = string.format("  f(x%d,yp)=f(%.4f,%.4f)=%.6f", i + 1, xnext, yp, fp)
        log[#log + 1] = string.format("  corrector y%d=y%d+(h/24)(9f(x%d,yp)+19f%d-5f%d+f%d)",
            i + 1, i, i + 1, i, i - 1, i - 2)
        log[#log + 1] = string.format("  y%d = %.6f   (x%d=%.4f)", i + 1, ynext, i + 1, xnext)

        xs[i + 1], ys[i + 1] = xnext, ynext
        fs[i + 1] = fnext
    end

    return xs, ys, log
end

local function compute()
    errorMsg, lines, scrollY = nil, nil, 0

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
    if n < 4 then
        errorMsg = "Se requieren al menos 4 pasos (n>=4): reduzca h o aumente xn"
        return
    end
    if n > MAX_N then
        errorMsg = "Demasiados pasos (maximo " .. MAX_N .. ")"
        return
    end

    local xs, ys, log, err = solve(expr, x0, y0, h, n)
    if not xs then
        errorMsg = err
        return
    end

    lines = log
    lines[#lines + 1] = ""
    lines[#lines + 1] = "TABLA DE RESULTADOS"
    lines[#lines + 1] = "  i        x            y"
    for i = 0, n do
        lines[#lines + 1] = string.format("  %-3d  %10s  %12s", i, fmt(xs[i]), fmt(ys[i]))
    end
end

function on.paint(gc)
    local w, hgt = platform.window:width(), platform.window:height()
    gc:setColorRGB(255, 255, 255)
    gc:fillRect(0, 0, w, hgt)
    gc:setColorRGB(0, 0, 0)

    gc:setFont("sansserif", "b", 10)
    gc:drawString("EDO: Adams-Bashforth-Moulton", 2, 2, "top")

    gc:setFont("sansserif", "r", 9)
    local y = 18
    local labelW = 42
    for i, f in ipairs(fields) do
        gc:drawString(f.label, 2, y, "top")
        local boxX, boxW = labelW, w - labelW - 4
        gc:setColorRGB(focus == i and 205 or 235, focus == i and 225 or 235, focus == i and 255 or 235)
        gc:fillRect(boxX, y, boxW, 14)
        gc:setColorRGB(0, 0, 0)
        gc:drawRect(boxX, y, boxW, 14)
        gc:drawString(f.value, boxX + 2, y, "top")
        y = y + 16
    end

    gc:setFont("sansserif", "r", 7)
    gc:drawString("Enter=calcular  Tab=cambiar campo  Arriba/Abajo=desplazar", 2, y, "top")
    y = y + 12

    gc:setFont("sansserif", "r", 8)
    if errorMsg then
        gc:setColorRGB(200, 0, 0)
        gc:drawString(errorMsg, 2, y, "top")
        gc:setColorRGB(0, 0, 0)
        return
    end

    if not lines then
        gc:drawString("Ingrese los datos y presione Enter", 2, y, "top")
        return
    end

    local visibleRows = math.floor((hgt - y) / ROW_H)
    local lastRow = math.min(#lines, scrollY + visibleRows)
    for r = scrollY + 1, lastRow do
        gc:drawString(lines[r], 2, y, "top")
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
    if not lines then return end
    local maxScroll = math.max(0, #lines - 1)
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
        if y >= yy and y <= yy + 14 then
            focus = i
            platform.window:invalidate()
            return
        end
        yy = yy + 16
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
