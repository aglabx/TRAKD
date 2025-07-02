# Makefile для проекта TRAKD
# Компилирует все C++ файлы в отдельные исполняемые файлы

# Компилятор и флаги
CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -pthread
LDFLAGS = -pthread

# Директории
SRCDIR = src
BINDIR = bin

# Исходные файлы
SOURCES = $(wildcard $(SRCDIR)/*.cpp)

# Исполняемые файлы (без расширения .cpp и с префиксом директории bin/)
EXECUTABLES = $(patsubst $(SRCDIR)/%.cpp,$(BINDIR)/%,$(SOURCES))

# Цель по умолчанию
all: $(BINDIR) $(EXECUTABLES)

# Создание директории bin
$(BINDIR):
	mkdir -p $(BINDIR)

# Правило для компиляции каждого .cpp файла в исполняемый файл
$(BINDIR)/%: $(SRCDIR)/%.cpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# Очистка скомпилированных файлов
clean:
	rm -rf $(BINDIR)

# Пересборка (очистка + сборка)
rebuild: clean all

# Установка зависимостей (если нужно)
install-deps:
	@echo "Для данного проекта специальные зависимости не требуются"
	@echo "Убедитесь, что у вас установлен g++ с поддержкой C++17"

# Запуск тестов (пример)
test: all
	@echo "Запуск базовых тестов..."
	@echo "Проверка компиляции завершена успешно"
	@ls -la $(BINDIR)/

# Показать справку
help:
	@echo "Доступные цели:"
	@echo "  all          - Скомпилировать все исполняемые файлы (по умолчанию)"
	@echo "  clean        - Удалить все скомпилированные файлы"
	@echo "  rebuild      - Очистить и пересобрать все"
	@echo "  install-deps - Показать информацию о зависимостях"
	@echo "  test         - Запустить простую проверку"
	@echo "  help         - Показать эту справку"
	@echo ""
	@echo "Исполняемые файлы будут созданы в директории bin/:"
	@echo "  bin/kmer_analyzer"
	@echo "  bin/locus_bed_generator" 
	@echo "  bin/LocusBedGeneratorDetailed"
	@echo "  bin/distance_analyzer_detailed"

# Дополнительные флаги для отладки
debug: CXXFLAGS += -g -DDEBUG
debug: all

# Объявление фиктивных целей
.PHONY: all clean rebuild install-deps test help debug
