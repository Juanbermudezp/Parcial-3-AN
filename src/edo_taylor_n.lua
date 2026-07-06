platform.apiLevel = '2.0'

-- Programa TI-Nspire: solucion numerica de EDO de primer orden
-- y' = f(x,y), y(x0) = y0, por el metodo de Taylor de orden n (n configurable).
-- Requiere una calculadora con CAS (TI-Nspire CX CAS / CX II CAS) porque
-- calcula las derivadas totales de f de forma simbolica.
-- Muestra el desarrollo paso a paso y al final la tabla de resultados.

local MAX_N = 100
local MAX_ORDER = 6
local ROW_H = 12
local SCROLL_STEP = 3

local fields = {
    { label = "f(x,y)=",  value = "x+y" },
    { label = "x0=",      value = "0" },
    { label = "y0=",      value = "1" },
    { label = "xn=",      value = "1" },
    { label = "h=",       value = "0.1" },
    { label = "orden n=", value = "2" },
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

-- Deriva simbolicamente 'expr' respecto de 'wrt' con el CAS y trae el
-- resultado a Lua como texto (string(...) fuerza el tipo string; sin CAS
-- esto falla porque el motor numerico no maneja expresiones simbolicas).
local function symbDeriv(expr, wrt)
    local res = math.eval("string(d(" .. expr .. "," .. wrt .. "))")
    if type(res) ~= "string" then
        return nil, "No se pudo derivar f(x,y) (se requiere calculadora con CAS)"
    end
    return res
end

local function symbSimplify(expr)
    local res = math.eval("string(" .. expr .. ")")
    if type(res) ~= "string" then return expr end
    return res
end

-- Construye T[1..n], las derivadas totales de y respecto de x sobre la
-- curva solucion:  T1 = f ;  Tk = d(T(k-1),x) + d(T(k-1),y) * f
local function buildTaylorCoeffs(expr, n)
    local T = { [1] = expr }
    local log = {}
    for k = 2, n do
        local dx, ex = symbDeriv(T[k - 1], "x")
        if not dx then return nil, nil, ex end
        local dy, ey = symbDeriv(T[k - 1], "y")
        if not dy then return nil, nil, ey end
        T[k] = symbSimplify("(" .. dx .. ")+(" .. dy .. ")*(" .. expr .. ")")
        log[#log + 1] = "  y^(" .. k .. ") = " .. T[k]
    end
    return T, log
end

local function solve(expr, x0, y0, h, n, order)
    local T, coeffLog, cerr = buildTaylorCoeffs(expr, order)
    if not T then return nil, nil, nil, cerr end

    local fact = { [0] = 1 }
    for k = 1, order do fact[k] = fact[k - 1] * k end

    local xs, ys = { [0] = x0 }, { [0] = y0 }
    local log = { "Derivadas totales usadas:", "  y' = " .. expr }
    for _, l in ipairs(coeffLog) do log[#log + 1] = l end
    log[#log + 1] = ""

    for i = 0, n - 1 do
        local terms, sum = {}, 0
        for k = 1, order do
            local val, verr = feval(T[k], xs[i], ys[i])
            if not val then return nil, nil, nil, verr end
            local term = (h ^ k / fact[k]) * val
            terms[k] = { val = val, term = term }
            sum = sum + term
        end
        local ynext = ys[i] + sum
        local xnext = xs[i] + h

        log[#log + 1] = string.format("Paso %d (x%d=%.4f,y%d=%.6f):", i + 1, i, xs[i], i, ys[i])
        for k = 1, order do
            log[#log + 1] = string.format("  y^(%d)=%.6f  termino=h^%d/%d!*y^(%d)=%.8f",
                k, terms[k].val, k, k, k, terms[k].term)
        end
        log[#log + 1] = string.format("  y%d = y%d + suma terminos = %.6f  (x%d=%.4f)",
            i + 1, i, ynext, i + 1, xnext)

        xs[i + 1], ys[i + 1] = xnext, ynext
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
    local order = toNumber(fields[6].value)

    if expr == "" or not x0 or not y0 or not xn or not h or not order then
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

    order = math.floor(order + 0.5)
    if order < 1 or order > MAX_ORDER then
        errorMsg = "orden n debe estar entre 1 y " .. MAX_ORDER
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

    local xs, ys, log, err = solve(expr, x0, y0, h, n, order)
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
    gc:drawString("EDO: Taylor de orden n", 2, 2, "top")

    gc:setFont("sansserif", "r", 9)
    local y = 18
    local labelW = 48
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
    y = y + 10
    gc:drawString("Nota: requiere calculadora con CAS", 2, y, "top")
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
