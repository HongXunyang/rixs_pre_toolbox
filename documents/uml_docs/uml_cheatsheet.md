# UML Quick Reference for RIXS Preparation Toolbox

This cheatsheet provides a quick reference for the UML notations most relevant to our project.

## Class Diagram Notation

```
┌───────────────────────┐
│      ClassName        │
├───────────────────────┤
│ - privateAttribute    │
│ + publicAttribute     │
├───────────────────────┤
│ + publicMethod()      │
│ - privateMethod()     │
└───────────────────────┘
```

### Relationships

```
ClassA ────────> ClassB     Association: ClassA uses ClassB
ClassA ◆────────> ClassB     Composition: ClassB is part of ClassA
ClassA ○────────> ClassB     Aggregation: ClassB is used by ClassA
ClassA ────────|> ClassB     Inheritance: ClassA inherits from ClassB
ClassA - - - -> ClassB     Dependency: ClassA depends on ClassB
```

## Sequence Diagram Notation

```
┌───┐          ┌───┐
│ A │          │ B │
└─┬─┘          └─┬─┘
  │              │
  │─────────────>│  Method call
  │              │
  │<- - - - - - -│  Return value
  │              │
  │──┐           │
  │  │ Self call │
  │<─┘           │
  │              │
```

## Component Diagram Notation

```
┌─────────────────┐       ┌─────────────────┐
│                 │       │                 │
│   Component A   │───────│   Component B   │
│                 │       │                 │
└─────────────────┘       └─────────────────┘
```

## PlantUML Quick Reference

### Class Diagram

```plantuml
@startuml
class MyClass {
  - privateAttribute
  + publicAttribute
  # protectedAttribute
  ~ packageAttribute
  + method(param: Type): ReturnType
}

interface MyInterface {
  + abstractMethod()
}

MyClass ..|> MyInterface: implements
@enduml
```

### Relationships

```plantuml
@startuml
ClassA --> ClassB: uses >
ClassA --* ClassB: contains
ClassA --o ClassB: aggregates
ClassA --|> ClassB: inherits from
ClassA ..> ClassB: depends on
@enduml
```

### Sequence Diagram

```plantuml
@startuml
actor User
participant "Class A" as A
participant "Class B" as B

User -> A: doSomething()
activate A
A -> B: process()
activate B
B --> A: result
deactivate B
A --> User: display result
deactivate A
@enduml
```

### Component Diagram

```plantuml
@startuml
package "Frontend" {
  [MainWindow]
  [Tabs]
}

package "Backend" {
  [Calculators]
  [Visualizers]
}

[MainWindow] --> [Tabs]
[Tabs] --> [Calculators]
[Tabs] --> [Visualizers]
@enduml
```

## Standard Tab Interface in Our Project

All new tabs should implement this interface:

```plantuml
@startuml
abstract class TabInterface {
  # main_window: MainWindow
  # layout: QGridLayout
  + __init__(main_window)
  + init_ui()
  + set_parameters(params: dict)
}

class YourNewTab {
  - calculator: YourCalculator
  - visualizer: YourVisualizer
  + __init__(main_window)
  + init_ui()
  + set_parameters(params: dict)
  + your_specific_methods()
}

TabInterface <|-- YourNewTab
@enduml
```

## Standard Calculator Interface

All calculators should implement this interface:

```plantuml
@startuml
class YourCalculator {
  - params: dict
  - initialized: bool
  + initialize(params: dict): bool
  + is_initialized(): bool
  + your_calculation_methods(): dict
}
@enduml
```

## Key Points to Remember

1. **Return Dictionaries**: All calculation methods should return dictionaries with at least:
   ```python
   {
       'success': bool,
       'error': str,  # Only if success is False
       # Other result data...
   }
   ```

2. **Standard Parameter Format**: Use consistent parameter structure:
   ```python
   params = {
       'lattice': {
           'a': float,
           'b': float,
           'c': float,
           'alpha': float,
           'beta': float,
           'gamma': float
       },
       'energy': float,
       # Other parameters...
   }
   ```

3. **Visualization Methods**: Always include options for axis and clearing:
   ```python
   def visualize_something(self, data, ax=None, clear=True):
       # Implementation
   ``` 