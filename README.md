# foamMultiRunCHT

# foamMultiRunCHT

## English

### Description

`foamMultiRunCHT` is a modified version of the classic `foamMultiRun` solver, developed for **OpenFOAM.org v14**.

The solver is designed for conjugate heat transfer (CHT) simulations involving two computational regions: `multicomponentFluid` (fluid region) and `solid` (solid region).

Its main feature is the ability to perform calculations in the fluid and solid regions using different time-stepping strategies and time-step sizes.

This allows the fluid and solid regions to be solved with different temporal resolutions, providing greater flexibility and computational efficiency for coupled simulations.

### Requirements and limitations

- OpenFOAM.org v14.
- The time-derivative numerical schemes (`ddtSchemes`) must use the **Euler** scheme in `fvSchemes`.

### Configuration

To configure the time-stepping parameters, add the following `regionTime` dictionary to `system/controlDict`:

```foam
regionTime
{
    fluid
    {
        mode            quasiSteady;
        deltaT          1e-6;
        maxIter         300;
        heatFluxTol     1e-3;
        heatFlux
        {
            region      fluid;
            patch       fluid_to_solid;
        }
    }

    solid
    {
        mode            transient;
        deltaT          1e-4;
        interval        1e-3;
    }
}
```

---

## Русский

### Описание

`foamMultiRunCHT` — модификация классического солвера `foamMultiRun`, разработанная для **OpenFOAM.org v14**.

Солвер предназначен для моделирования сопряжённого теплообмена (Conjugate Heat Transfer, CHT) между двумя вычислительными регионами: `multicomponentFluid` (жидкостная область) и `solid` (твёрдое тело).

Основная особенность солвера — возможность выполнять расчёты в жидкостной и твёрдой областях с использованием различных стратегий временного интегрирования и разных шагов по времени.

Это позволяет задавать различную временную дискретизацию для каждого региона, обеспечивая большую гибкость и эффективность вычислений при моделировании сопряжённых задач.

### Требования и ограничения

- OpenFOAM.org v14.
- Численные схемы для производной по времени (`ddtSchemes`) в файле `fvSchemes` должны использовать схему **Euler**.

### Настройка

Для настройки параметров расчёта необходимо добавить следующий словарь `regionTime` в файл `system/controlDict`:

```foam
regionTime
{
    fluid
    {
        mode            quasiSteady;
        deltaT          1e-6;
        maxIter         300;
        heatFluxTol     1e-3;
        heatFlux
        {
            region      fluid;
            patch       fluid_to_solid;
        }
    }

    solid
    {
        mode            transient;
        deltaT          1e-4;
        interval        1e-3;
    }
}
```
